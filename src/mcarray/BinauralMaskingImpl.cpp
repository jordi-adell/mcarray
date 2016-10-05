/*
* BinauralMaskingImpl.cpp
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

#include <mcarray/BinauralMaskingImpl.h>
#include <mcarray/mcalogger.h>
#include <mcarray/mcarray_exception.h>
#include <mcarray/mcadefs.h>
#include <mcarray/microhponeArrayHelpers.h>


#include <wipp/wipputils.h>
#include <wipp/wippstats.h>
#include <wipp/wippsignal.h>

#include <sstream>

#include <math.h>

namespace mca {

// nBins + 1 is to store the residual of the FB.
BinauralMaskingImpl::BinauralMaskingImpl(int samplerate, double microDistance, float lowFreq, float highFreq, MaskingMethod mmethod) :
  ShortTimeProcess(calculateWindowSizeFromSampleRate(samplerate, _frameRate),
		   calculateWindowSizeFromSampleRate(samplerate, _frameRate)*(_nBins+1),
		   2),
  _fftOrder(calculateOrderFromSampleRate(samplerate, _frameRate)),
  _sampleRate(samplerate),
  _microDistance(microDistance),
  _mmethod(mmethod),
  _minFreq(lowFreq),
  _maxFreq(highFreq),
  _nchannels(getNumberOfChannels()),
  _windowSize(getFrameSize())
{
  init();
  TRACE_STREAM("Binaural Localisation parameters: " << std::endl
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


BinauralMaskingImpl::~BinauralMaskingImpl()
{

}

void BinauralMaskingImpl::init()
{
  if (_nchannels != 2)
  {
    throw(MCArrayException("Sound localisation is only working for 2 channels by now."));
  }
  _filterBankRight.reset(new dsp::FilterBankFFTWMelScale(_fftOrder, _nBins, _sampleRate, _minFreq, _maxFreq));
  _filterBankLeft.reset(new dsp::FilterBankFFTWMelScale(_fftOrder, _nBins, _sampleRate, _minFreq, _maxFreq));
  _shortTimePower.reset(new BaseType[_nBins]);
  wipp::setZeros(_shortTimePower.get(), _nBins);
  calculateThresholds();
}


void BinauralMaskingImpl::frameAnalysis(BaseType *inFrame, BaseType *analysis, int frameLength, int analysisLength, int channel)
{
  dsp::FilterBank *filter;
  if (channel == 0)
    filter = _filterBankLeft.get();
  else
    filter = _filterBankRight.get();
  filter->filterBuffer(inFrame, &analysis[_nBins*frameLength] , analysis, frameLength, analysisLength);
}


void BinauralMaskingImpl::processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
						 std::vector<double*> &dataChannels, int dataLength)
{
  BaseType *left =  analysisFrames[0];
  BaseType *right = analysisFrames[1];
  int nspatialMaskedBins = 0;
  int ntempMaskedBins = 0;
  int nenhancedBins= 0;
  std::ostringstream oss;
  //  oss << std::setprecision(4);

  oss << "[";
  for (int bin = 0; bin < _nBins; ++bin) // Residual should not be processed (i<_nBins) instead of (i<=_nBins)
  {
    BaseType *binleft  = &left[bin*_windowSize];
    BaseType *binright = &right[bin*_windowSize];

    // temporal masking needs to be executed always, because it has memeory
    bool tempMask = temportalMasking(binleft, binright, _windowSize, bin);
    bool spatMask = spatialMasking(binleft, binright, _windowSize, bin);

    if (spatMask)
    {
      maskFrame(binleft, _windowSize, _spatialMaskingFactor, bin);
      maskFrame(binright, _windowSize, _spatialMaskingFactor, bin);
      oss << _filterBankLeft.get()->getBinCenterFrequency(bin)*_sampleRate << " ";
      ++nspatialMaskedBins;

    }
    else if (tempMask)
    {
      maskFrame(binleft, _windowSize, _temporalMaskingFactor, bin);
      maskFrame(binright, _windowSize, _temporalMaskingFactor, bin);
      oss << _filterBankLeft.get()->getBinCenterFrequency(bin)*_sampleRate <<  " ";
      ++ntempMaskedBins;
    }
    else
    {
      enhanceFrame(binleft, _windowSize);
      enhanceFrame(binright, _windowSize);
      ++nenhancedBins;
    }
  }
  //          TRACE_STREAM("SM: " << nspatialMaskedBins << " bins, "
  //                           << "TM: " << ntempMaskedBins << " bins, "
  //                           << "EB: " << nenhancedBins << " bins, "
  //                           << "total: " << _nBins << " bins");
  oss << "]";
  //          TRACE_STREAM(oss.str());
}

void BinauralMaskingImpl::frameSynthesis(BaseType *outFrame, BaseType *analysis, int frameLength, int analysisLength, int)
{
  wipp::setZeros(outFrame, frameLength);
  for (int bin = 0, offset = 0;
       (bin <= _nBins) && (offset < analysisLength - frameLength);
       ++bin, offset += frameLength)
  {
    wipp::add(&analysis[offset], outFrame, frameLength);
  }
}


void BinauralMaskingImpl::zeroFrame(BaseType *frame, int length)
{
  wipp::divC(1000, frame, length); // 1/1000 --> -60dB
}

void BinauralMaskingImpl::maskFrameByScaling(BaseType *frame, int length, int bin) // RELATIVE
{

  //
  //                ro * (1/N) *  sum_n(frame[n]²)
  // factor = sqrt( ------------------------------ )
  //                           Q[m-1]
  //
  // ro = scalingFactor
  // Q[m-1] - low-pass filtered power in the previous frame
  //

  BaseType aux[length];
  BaseType factor;
  wipp::sqr(frame, aux, length);  // frame[n]^2
  wipp::mean(aux, length, &factor);   // (1/N)*sum_n
  factor *= _scalingFactor;   // =* ro
  factor /= _shortTimePower[bin];   //  =* (1/Q[m-1])
  factor = sqrt(factor);   // sqrt
  wipp::multC(factor, frame, length);
}

void BinauralMaskingImpl::maskFrameByFactor(BaseType *frame, int length, float factor)
{
  wipp::divC(factor, frame, length);
}

void BinauralMaskingImpl::maskFrame(BaseType *frame, int length,float factor, int bin)
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
  }
}

void BinauralMaskingImpl::enhanceFrame(BaseType *frame, int length)
{
  wipp::multC(_enhanceFactor, frame, length);
}

double BinauralMaskingImpl::localise(BaseType *left, BaseType *right, int length)
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

      //      ippsCrossCorr_64f(&left[offset], _windowSize, &right[offset], _windowSize, crossCorr, _windowSize, _windowSize/2);
      wipp::cross_corr(&left[offset], _windowSize, &right[offset], _windowSize, crossCorr, _windowSize, _windowSize/2);
      wipp::maxidx(crossCorr, _windowSize, &max, &idx);
      angle = (static_cast<double>(idx) - _windowSize/2)*2*M_PI/_windowSize;
    }
  }
  return angle;
}


void BinauralMaskingImpl::calculateThresholds()
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
    double wfreq =  _filterBankLeft.get()->getBinCenterFrequency(bin)*_sampleRate*2*M_PI;
    _thresholds[bin] = cos(wfreq*_microDistance*sin(_phi)/getSpeedOfSound()) * 0.9;
    TRACE_STREAM("BIN: " << bin << " F: " << wfreq/(2*M_PI) << " TH: " << _thresholds[bin]);
  }

}

// returns true if the signal has to be masked (removed)
inline bool BinauralMaskingImpl::spatialMasking(BaseType *left, BaseType *right, int length, int bin)
{
  double ncorr = normaliseCorrelation(left, right, length);
  //          TRACE_STREAM("bin: " << bin << " f: " << _filterBankLeft.get()->getBinCenterFrequency(bin)*_sampleRate << " -->  corr: " << ncorr << " vs. " <<  _thresholds.at(bin));
  bool masking =  ( ncorr < _thresholds.at(bin));
  return masking;
}

double BinauralMaskingImpl::normaliseCorrelation(BaseType *left, BaseType *right, int length)
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

  wipp::sqr(right, vaux, length);
  wipp::mean(vaux, length, &aux);
  denom *= sqrt(aux);

  if (denom == 0)
    return 1;
  else
    return numer/denom;
}



inline bool BinauralMaskingImpl::temportalMasking(BaseType *left, BaseType *right, int length, int bin)
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
  return (power < _shortTimePower[bin]);
}

inline double BinauralMaskingImpl::getFramePower(BaseType *left, BaseType *right, int length)
{

  //
  //  (1/N) * sum_n( ((left[n]+right[n])/2)² )
  //

  BaseType vaux[length];
  BaseType power;
  wipp::copyBuffer(left, vaux, length);
  wipp::add(right, vaux, length);
  wipp::divC(2, vaux, length);
  wipp::sqr(vaux, length);
  wipp::mean(vaux, length, &power);
  return power;
}


}
