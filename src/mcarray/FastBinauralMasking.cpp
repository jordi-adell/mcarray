/*
* FastBinauralMasking.cpp
* Copyright 2016 (c) Jordi Adell
* Created on: 2015
* 	Author: Jordi Adell - adellj@gmail.com
*
* This file is part of MCARRAY
*
* MCARRAY is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* MCARRAY is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with MCARRAY.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <mcarray/FastBinauralMasking.h>
#include <mcarray/mcalogger.h>
#include <mcarray/mcarray_exception.h>

#include <dspone/rt/ShortTimeFourierTransform.h>
#include <dspone/filter/FilterBankMelScale.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>

#include <sstream>
#include <iomanip>
#include <math.h>


#define BOTH     BinauralMasking::BOTH
#define NOTHING  BinauralMasking::NOTHING
#define SPATIAL  BinauralMasking::SPATIAL
#define FULL     BinauralMasking::FULL
#define RELATIVE BinauralMasking::RELATIVE
#define FACTOR   BinauralMasking::FACTOR
#define NOISY    BinauralMasking::NOISY


namespace mca {


// nBins + 1 is to store the residual of the FB.
FastBinauralMasking::FastBinauralMasking(int samplerate,
					 double microDistance,
					 float lowFreq,
					 float highFreq,
					 MaskingMethod mmethod,
					 MaskingAlg algorithm) :
    dsp::STFT(2 ,calculateOrderFromSampleRate(samplerate, _frameRate)),
    _sampleRate(samplerate),
    _microDistance(microDistance),
    _mmethod(mmethod),
    _algorithm(algorithm),
    _minFreq(lowFreq),
    _maxFreq(highFreq),
    _nchannels(getNumberOfChannels()),
    _fftOrder(calculateOrderFromSampleRate(samplerate, _frameRate)),
    _oneSidedFFTLength(getOneSidedFFTLength()),
    _windowSize(getWindowSize()),
    _firstCall(0)
{
    init();
    TRACE_STREAM("FAST Binaural Localisation parameters: " << std::endl
		 << "order: "  << _fftOrder
		 << ", rate: " << _sampleRate
		 << ", micro distance: " << _microDistance
		 << ", method: " << _mmethod << std::endl
		 << "temporal factor: " << _temporalMaskingFactor
		 << ", spatial factor: " << _spatialMaskingFactor
		 << ", scaling factor: " << _scalingFactor
		 << ", enhance factor: " << _enhanceFactor
		 << ", Bandwidth: (" << _minFreq << ", " << _maxFreq << ")"
		 );
}


void FastBinauralMasking::init()
{

    if (_nchannels != 2)
    {
	throw(MCArrayException("Binaural masking is only working for 2 channels."));
    }

	      //    _coeficientsLength = _nBins * ( (1 << (_fftOrder-1)) + 1 );
    _coeficientsLength = getAnalysisLength() * _nBins;
    _filterBank.reset(new dsp::FilterBankFFTWMelScale(_fftOrder, _nBins, _sampleRate, _minFreq, _maxFreq));
    _filterCoeficients.reset(new BaseTypeC[_coeficientsLength]);
    BaseType coefs[_coeficientsLength];
    int gotNCoefs = _filterBank->getFiltersCoeficients(coefs, _coeficientsLength);
    if (gotNCoefs != _coeficientsLength)
    {
	std::ostringstream oss;
	oss << "Number of filter coeficients different from expected: " << gotNCoefs << " instead of " << _coeficientsLength;
	throw(MCArrayException(oss.str()));
    }
    DEBUG_STREAM("number of coefs: " << _coeficientsLength << " order: " << _fftOrder);
    DEBUG_STREAM("BW: " << _minFreq << " - " << _maxFreq);
    //    ippsRealToCplx_64f(coefs, NULL, _filterCoeficients.get(), _coeficientsLength);
    wipp::real2complex(coefs, NULL, reinterpret_cast<wipp::wipp_complex_t*>(_filterCoeficients.get()), _coeficientsLength);

    _shortTimePower.reset(new BaseType[_nBins]);
    _noiseEstimatePower.reset(new BaseType[_nBins]);
    _fftLeftFrame.reset(new BaseType[getAnalysisLength()]);
    _fftRightFrame.reset(new BaseType[getAnalysisLength()]);
    _outLeftFrame.reset(new BaseType[getAnalysisLength()]);
    _outRightFrame.reset(new BaseType[getAnalysisLength()]);

    wipp::setZeros(_shortTimePower.get()    , _nBins);
    wipp::setZeros(_noiseEstimatePower.get(), _nBins);

    //          _gcc.reset(new GeneralisedCrossCorrelation(_analysisLength ,GeneralisedCrossCorrelation::ONESIDEDFFT));

    calculateThresholds();
}


void FastBinauralMasking::processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
						 std::vector<double*> &dataChannels, int dataLength)
{

    if (_mmethod == NOTHING)
    {
	WARN_STREAM_ONCE("Fast Binaural masking is disabled.");
	return;
    }

    int nspatialMaskedBins = 0;
    int ntempMaskedBins = 0;
    int nenhancedBins= 0;
    std::ostringstream oss;
    oss << std::setprecision(4);

    wipp::setZeros(_outLeftFrame.get(), analysisLength);
    wipp::setZeros(_outRightFrame.get(), analysisLength);

    oss << "[";
    for (int bin = 0; bin < _nBins; ++bin) // Residual should not be processed (i<_nBins) instead of (i<=_nBins)
    {
      wipp::mult(reinterpret_cast<wipp::wipp_complex_t*>(analysisFrames[0]),
	  reinterpret_cast<wipp::wipp_complex_t*>(&_filterCoeficients.get()[bin*_oneSidedFFTLength]),
	  reinterpret_cast<wipp::wipp_complex_t*>(_fftLeftFrame.get()),  _oneSidedFFTLength);
      wipp::mult(reinterpret_cast<wipp::wipp_complex_t*>(analysisFrames[1]),
	  reinterpret_cast<wipp::wipp_complex_t*>(&_filterCoeficients.get()[bin*_oneSidedFFTLength]),
	  reinterpret_cast<wipp::wipp_complex_t*>(_fftRightFrame.get()), _oneSidedFFTLength);

	// temporal masking needs to be executed always, because it has memory
	bool tempMask, spatMask=false;

	tempMask = temportalMasking(_fftLeftFrame.get(),  _fftRightFrame.get(), _windowSize, bin);
	if (_algorithm == BOTH || _algorithm == SPATIAL)
	{
	    spatMask = spatialMasking(_fftLeftFrame.get(),  _fftRightFrame.get(), _windowSize, bin);
	    if (_algorithm == SPATIAL)
	    {
		tempMask = false;
	    }
	}

	if (spatMask)
	{
	    maskFrame(_fftLeftFrame.get(),  _windowSize, _spatialMaskingFactor, bin);
	    maskFrame(_fftRightFrame.get(), _windowSize, _spatialMaskingFactor, bin);
	    oss << _filterBank.get()->getBinCenterFrequency(bin)*_sampleRate << " ";
	    ++nspatialMaskedBins;

	}
	else if (tempMask)
	{
	    maskFrame(_fftLeftFrame.get(),  _windowSize, _temporalMaskingFactor, bin);
	    maskFrame(_fftRightFrame.get(), _windowSize, _temporalMaskingFactor, bin);
	    oss << _filterBank.get()->getBinCenterFrequency(bin)*_sampleRate <<  " ";
	    ++ntempMaskedBins;
	}
	else
	{
	    enhanceFrame(_fftLeftFrame.get(),  _windowSize);
	    enhanceFrame(_fftRightFrame.get(), _windowSize);
	    ++nenhancedBins;
	}
	wipp::add(_fftLeftFrame.get(), _outLeftFrame.get(), analysisLength);
	wipp::add(_fftRightFrame.get(), _outRightFrame.get(), analysisLength);
    }

    ++_firstCall;
    if (_firstCall < 2)
    {
	wipp::copyBuffer(_shortTimePower.get(), _noiseEstimatePower.get(), _nBins);
    }

    wipp::copyBuffer(_outLeftFrame.get(),  analysisFrames[0], analysisLength);
    wipp::copyBuffer(_outRightFrame.get(), analysisFrames[1], analysisLength);



    TRACE_STREAM("SM: " << nspatialMaskedBins << " bins, "
		 << "TM: " << ntempMaskedBins << " bins, "
		 << "EB: " << nenhancedBins << " bins, "
		 << "total: " << _nBins << " bins");
    oss << "]";
    TRACE_STREAM(oss.str());
}



void FastBinauralMasking::zeroFrame(BaseType *frame, int length)
{
    wipp::divC(1000, frame, length); // 1/1000 --> -60dB
}


void FastBinauralMasking::noisyFrame(BaseType *frame, int length, int bin)
{
    BaseType powerBin = getPower(frame,length);
    BaseType noisePower = _noiseEstimatePower[bin];
    BaseType factor = 1;
    if (powerBin > 0)
	factor = noisePower/powerBin;
    //            BaseType mean=100;
    SignalPtr constant;
    constant.reset(new BaseType[_oneSidedFFTLength]);

    if (_firstCall < 2)
	return;

    //            ippsMean_64f(_noiseEstimatePower.get(),_nBins, &mean);
    //            ippsSet_64f(mean, constant.get(), _oneSidedFFTLength);
    //            ippsRealToCplx_64f(constant.get(), NULL, reinterpret_cast<BaseTypeC*>(frame), _oneSidedFFTLength);
    //            ippsMul_64fc_I(&_filterCoeficients.get()[bin*_oneSidedFFTLength],
    //                         reinterpret_cast<BaseTypeC*>(frame), _oneSidedFFTLength);

    //            ippsSet_64f(mean, frame, length*_maxFreq/_sampleRate);


    wipp::multC(factor, frame, length);

    //            ippsDiv_64f(&powerBin, &meanPower, &SNR, 1);

    //            factor=5.623/SNR;             //5.623 --> 15 dB

    //            wipp::multC(factor, frame, length);

    //            if(factor>0)
    //                wipp::divC(factor, frame, length);
    //            else
    //                wipp::multC(-factor, frame, length);

}



void FastBinauralMasking::maskFrameByScaling(BaseType *frame, int length, int bin) // RELATIVE
{

    //
    //                ro * (1/N) *  sum_n(frame[n]²)
    // factor = sqrt( ------------------------------ )
    //                           Q[m-1]
    //
    // ro = scalingFactor
    // Q[m-1] - low-pass filtered power in the previous frame
    //

    BaseType aux[_oneSidedFFTLength];
    BaseType factor;

    //I use magnitude because power is calculated in the Freq domain.
    wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(frame), aux, _oneSidedFFTLength); //  |·|
    wipp::sqr(aux, _oneSidedFFTLength); // frame[n]^2
    wipp::mean(aux, _oneSidedFFTLength, &factor);   // (1/N)*sum_n

    factor *= _scalingFactor;   // =* ro
    if (_shortTimePower[bin] < 1e-10)
	factor = _scalingFactor;
    else
	factor /= _shortTimePower[bin];   //  =/ Q[m-1]
    factor = sqrt(factor);   // sqrt
    wipp::multC(factor, frame, length);
}

void FastBinauralMasking::maskFrameByFactor(BaseType *frame, int length, float factor)
{
    wipp::divC(factor, frame, length);
}

void FastBinauralMasking::maskFrame(BaseType *frame, int length,float factor, int bin)
{
    switch(_mmethod)
    {
	case FULL:
	    zeroFrame(frame, length);
	break;
	case RELATIVE:
	    maskFrameByScaling(frame, length, bin);
	break;
	case FACTOR:
	    maskFrameByFactor(frame, length, factor);
	break;
	case NOISY:
	    noisyFrame(frame, length, bin);
	break;
	case NOTHING:
	break;
    }
}

void FastBinauralMasking::enhanceFrame(BaseType *frame, int length)
{
    wipp::multC(_enhanceFactor, frame, length);
}

double FastBinauralMasking::localise(BaseType *left, BaseType *right, int length)
{
    double angle = 0;
    for (int channel = 0; channel < 2; ++channel)
    {
	for (int bin = 0, offset = 0;
	     offset < length &&  bin < _nBins;
	     ++bin, offset+=_windowSize)
	{
	    BaseType crossCorr[_windowSize];
	    BaseType max;
	    size_t idx;

	    wipp::cross_corr(&left[offset], _windowSize, &right[offset], _windowSize, crossCorr, _windowSize, _windowSize/2);
	    wipp::maxidx(crossCorr, _windowSize, &max, &idx);
	    angle = (static_cast<double>(idx) - _windowSize/2)*2*M_PI/_windowSize;
	}
    }
    return angle;
}


void FastBinauralMasking::calculateThresholds()
{

    //
    // cos(w_0*tau) * cos(w_0*(d/c)*sin(_phi))
    //
    // tau - delay between microphones
    // d*sin(_phi) delay distance (in metres)
    // d*sin(_phi)/c delay in seconds
    //
    // w_0 = 2*pi*f_0
    // (where f is the analog frequency at the center of the bin)
    //
    // d/c - where d = distance between micros and c = speed of sound
    //

    _thresholds.resize(_nBins, 0);
    for (int bin = 0; bin<_nBins; ++bin)
    {
	double wfreq =  _filterBank.get()->getBinCenterFrequency(bin)*_sampleRate*2*M_PI;
	_thresholds[bin] = cos(wfreq*_microDistance*sin(_phi)/getSpeedOfSound());
	TRACE_STREAM("BIN: " << bin << " F: " << wfreq/(2*M_PI) << " TH: " << _thresholds[bin]);
    }

}

// returns true if the signal has to be masked (removed)
inline bool FastBinauralMasking::spatialMasking(BaseType *left, BaseType *right, int length, int bin)
{
    //          double ncorr = normaliseFFTCorrelation(left, right, length);
    double ncorr = generalisedCrossCorrelation(left, right, length);
    TRACE_STREAM("bin: " << bin << " f: " << _filterBank.get()->getBinCenterFrequency(bin)*_sampleRate << " -->  corr: " << ncorr << " vs. " <<  _thresholds.at(bin));
    bool masking =  ( ncorr < _thresholds.at(bin));
    return masking;
}

double FastBinauralMasking::normaliseCorrelation(BaseType *left, BaseType *right, int length)
{

    //                            (1/N)*sum_n(left[n]*right[n])
    //        --------------------------------------------------------------------------
    //         sqrt((1/N)*sum_n(left[n]*left[n]))*sqrt((1/N)*sum_n(right[n]*right[n]))

    BaseType numer;
    BaseType denom;
    BaseType aux;
    BaseType vaux[length];

    wipp::mult(left, right, vaux, length);
    wipp::mean(vaux, length, &numer);

    wipp::sqr(left, vaux, length);
    wipp::mean(vaux, length, &aux);
    denom = sqrt(aux);

    wipp::sqrt(right, vaux, length);
    wipp::mean(vaux, length, &aux);
    denom *= sqrt(aux);

    if (denom == 0)
	return 1;
    else
	return numer/denom;
}




double FastBinauralMasking::normaliseFFTCorrelation(BaseType *left, BaseType *right, int length)
{

    // This function calculates normalised correlation in the frequency domain,
    // Using the cross-correlation definition in freq domain and the Parseval's theorem.
    //                    sum_w [ W(w)X_i(w)X*_j(w) ]
    //        GCC  =     ------------------------------
    //                                 N

    //                                      N^2
    //        W(w) =  sqrt(-----------------------------------------)  ==>  PHAT (more or less)
    //                       sum_w(|X_i(w)|^2)*sum_w(|X_j(w)|^2)

    if (length != 2*(_oneSidedFFTLength-1))
    {
	std::ostringstream oss;
	oss << "The length of the frame to process do not match the one-sised FFT length: " << length
	    << "  vs. " << 2*_oneSidedFFTLength << " ("<< _fftOrder << ")";
	throw(MCArrayException(oss.str()));
    }

    BaseTypeC vauxC[_oneSidedFFTLength];
    BaseType vaux[_oneSidedFFTLength];
    wipp::wipp_complex_t meanC;
    BaseType numer, denom, mean;

    // Cross-Correlation
    wipp::conj(reinterpret_cast<wipp::wipp_complex_t*>(left), reinterpret_cast<wipp::wipp_complex_t*>(vauxC), _oneSidedFFTLength);
    wipp::mult(reinterpret_cast<wipp::wipp_complex_t*>(right), reinterpret_cast<wipp::wipp_complex_t*>(vauxC), _oneSidedFFTLength);
    wipp::mean(reinterpret_cast<wipp::wipp_complex_t*>(vauxC), _oneSidedFFTLength, &meanC);
    wipp::real(&meanC, &numer, 1);

    if (numer == 0)
	return 0;

    // Energy from Parseval's theorem
    wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(left), vaux, _oneSidedFFTLength);
    wipp::sqr(vaux, _oneSidedFFTLength);
    wipp::mean(vaux, _oneSidedFFTLength, &denom);

    wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(right), vaux, _oneSidedFFTLength);
    wipp::sqr(vaux, _oneSidedFFTLength);
    wipp::mean(vaux, _oneSidedFFTLength, &mean);

    denom = denom * mean;
    denom = sqrt(denom);
    if (denom == 0)
	return 1;
    else
	return numer/denom;
}


double FastBinauralMasking::generalisedCrossCorrelation(BaseType *left, BaseType *right, int length)
{
    //          BaseTypeC corr =  _gcc->calculateCorrelation(reinterpret_cast<BaseTypeC*>(left),
    //                                                       reinterpret_cast<BaseTypeC*>(right),
    //                                                       _fft.getOneSidedFFTLength(), 0,
    //                                                       GeneralisedCrossCorrelation::ONESIDEDFFT);
    //          return corr.re;

    double corr = normaliseFFTCorrelation(left, right, length);
    return corr;

}


inline bool FastBinauralMasking::temportalMasking(BaseType *left, BaseType *right, int length, int bin)
{

    //
    // first order IIR low pass filter:
    // Q[m] = labda*Q[m-1] + (1-lambda)*P[m]
    //
    // lambda = forgetingFactor
    // Q[m-1] - previous short-time power
    // P[m] - current frame power
    //

    BaseType power = getFramePower(left, right, length);
    _shortTimePower[bin] = _shortTimePower[bin]*_forgetingFactor + (1 - _forgetingFactor)*power;
    TRACE_STREAM("Q[m]=" << _shortTimePower[bin] << ", " << power);
    return (power < _rejectTemporalFactor*_shortTimePower[bin]);
}


inline double FastBinauralMasking::getFramePower(BaseType *left, BaseType *right, int length)
{

    //
    //  (1/N) * sum_n( ((left[n]+right[n])/2)² )
    //

    //          if (length != 2*(_oneSidedFFTLength-1))
    //            throw(SpeechException("Mismatch in vector lengths"));

    BaseType halfRight[length];
    BaseType halfLeft[length];
    BaseType mixed[length];

    wipp::copyBuffer(left, halfRight, length);
    wipp::copyBuffer(right, halfLeft, length);
    wipp::divC(2, halfLeft, length);
    wipp::divC(2, halfRight, length);
    wipp::add(halfRight, halfLeft, mixed, length);

    return getPower(mixed, length);
}



inline double FastBinauralMasking::getPower(BaseType *frame, int length)
{
    //
    //  (1/N) * sum_n( (frame[n])² )
    //

    int complexLength = length/2;

    BaseType magnitude[complexLength];
    BaseType power = 0;

    wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(frame), magnitude, complexLength);
    wipp::sqr(magnitude, complexLength);
    wipp::mean(magnitude, complexLength, &power);
    power = sqrt(power);
    return power;

}

}


