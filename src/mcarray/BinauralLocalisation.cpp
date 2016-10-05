/*
* BinauralLocalisation.cpp
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
#include <mcarray/BinauralLocalisation.h>
#include <mcarray/SoundLocalisationCallback.h>
#include "SoundLocalisationParticleFilter.h"

#include <dspone/algorithm/signalPower.h>

#include <wipp/wippsignal.h>
#include <wipp/wipputils.h>
#include <wipp/wippstats.h>

#include <math.h>
#include <sstream>


// This is because the use of a particle filter in localisation
// is not properly tested.
#define USE_PARTICLE_FILTER
//#undef  USE_PARTICLE_FILTER


namespace mca
{


//----------------- Temporal GCC Binaural Localisation ---------------------------------------------------------------

TemporalGCCBinauralLocalisation::TemporalGCCBinauralLocalisation(int sampleRate, ArrayDescription microphonePositions) :
    SoundLocalisationImpl(microphonePositions),
    ShortTimeAnalysis(static_cast<int>(2*(_frameRate*sampleRate)), static_cast<int>(2*(_frameRate*sampleRate))  , 2),
    _microphoneDistance(microphonePositions.distance(0,1)),
    _sampleRate(sampleRate),
    _ndelays(static_cast<int>(_microphoneDistance*_sampleRate/getSpeedOfSound())),
    _windowSize(getWindowSize()),
    _analysisLength(getAnalysisLength()),
    _signalPower(0)
{
    if (_ndelays < 2)
    {
	std::ostringstream oss;
	oss << "Sample frequency has to be a least " << getSpeedOfSound()*2/_microphoneDistance;
    }
    if (_microphonePositions.size() != 2)
    {
	WARN_STREAM("The number of microphones in ArrayDescription is different from 2.");
    }

    _currentDOA.reset(new BaseType[1]);
    _prob.reset(new BaseType[1]);
    _currentDOA[0] = 0;
    _prob[0] = -1;

    _leftDelayed.reset(new BaseType[_windowSize]);
    _rightDelayed.reset(new BaseType[_windowSize]);
    _crosCorr.reset(new BaseType[2*_ndelays+1]);
    _index.reset(new BaseType[_ndelays]);
    _mixedChannel.reset(new BaseType[_analysisLength]);
    _triangle.reset(new BaseType[_ndelays]);
    double phase = 3*M_PI_2;
    wipp::triangle(_triangle.get(), _ndelays, 1/(2*static_cast<double>(_ndelays)), phase, 0, 0);
    wipp::multC(0.1, _triangle.get(), _ndelays);
    //    ippsTriangle_Direct_64f(_triangle.get(), _ndelays, 0.1, 1/(2*static_cast<double>(_ndelays)), 0, &phase);
}
void TemporalGCCBinauralLocalisation::frameAnalysis(BaseType *inFrame, BaseType *analysis, int frameLength, int analysisLength, int)
{
    wipp::setZeros(analysis, analysisLength);
    //    wipp::copyBuffer(inFrame, analysis, std::min(frameLength, analysisLength));
    //    wipp::mult(_iwindow.get(), analysis, std::min(frameLength, analysisLength)); // To unwindow the frame.
    unwindowFrame(inFrame, analysis, std::min(frameLength, analysisLength));
}

//  log(_index[maxindex] / sum_i abs(_index[i]))
void TemporalGCCBinauralLocalisation::estimateDOALogLikelihood(int maxindex, BaseType &prob)
{
    BaseType aux[_ndelays];
    BaseType probFloor = 0;
    prob = log (0.0000000001);

    if (maxindex < _ndelays)
    {
      wipp::abs(_index.get(), aux, _ndelays);
      //	ippsAbs_64f(_index.get(), aux, _ndelays);
      wipp::sum(aux, _ndelays, &probFloor);
      probFloor += 0.000001;
      if (_index[maxindex] != 0)
      {
	prob = log(_index[maxindex] / probFloor);
      }
    }
    else
    {
      prob = 0;
    }
}

void TemporalGCCBinauralLocalisation::processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
							     std::vector<double*> &dataChannels, int dataLength)
{
    BaseType DOA = 0;
    BaseType power = 0;
    BaseType max = 0;

    // Getting signal power, and setting as floor if the first signal frame is processed.
    if (!_noiseEstimated)
    {
	power = setPowerFloor(analysisFrames, analysisLength, 2, _sampleRate);
    }
    else
    {
	power = dsp::SignalPower::logPower(analysisFrames, analysisLength);
    }

    DEBUG_STREAM("power " << power);

    // Only if singal power exceeds the floor power the DOA is obtained
    if (power > _powerFloor)
    {
	for (int nleftDelay = 0, nrightDelay= _ndelays;
	     nleftDelay < _ndelays && nrightDelay > 0;
	     ++nleftDelay, --nrightDelay)
	{
	    delay(analysisFrames[0], analysisLength, _leftDelayed.get(),  _windowSize, nleftDelay);
	    delay(analysisFrames[1], analysisLength, _rightDelayed.get(), _windowSize, nrightDelay);

	    crossCorrelation(_leftDelayed.get(), _windowSize, _rightDelayed.get(), _windowSize, _crosCorr.get(), 2*_ndelays+1);
	    _index[nleftDelay] = integration(_crosCorr.get(), 2*_ndelays+1);
	    // TODO this approach is much faster and I guess it is more accurate also, but needs testing.
	    //              _index[nleftDelay] = correlation(_leftDelayed.get(), _rightDelayed.get(), _windowSize);
	}

	wipp::add(_triangle.get(), _index.get(), _ndelays);
	maxminNormalisation(_index.get(), _ndelays);
	size_t maxIndex = 0;
	wipp::maxidx(_index.get(), _ndelays, &max, &maxIndex);
	DOA = samples2Degrees(maxIndex, _ndelays) - 90;
	estimateDOALogLikelihood(maxIndex, _prob[0]);
	_currentDOA[0] = _currentDOA[0]*_doaMemoryFactor + (1 - _doaMemoryFactor)*DOA;
    }
    else
    {
	// if the signal power is not high enough DOA is cosidered 0º (strict front)
	_currentDOA[0] = _currentDOA[0]*_doaMemoryFactor + (1 - _doaMemoryFactorSilence)*0;
	_prob[0] = -100000;
    }

    // The current DOA is updated considered past values with a memory factor.


    // Only if singal power exceeds the floor power the DOA is updated and
    // only in the same situation
    if (power > _powerFloor)
	_ptrCallback->setDOA(_currentDOA, _prob, power,1);

}

// outFrame[m] = inFrame[m-delay]
void TemporalGCCBinauralLocalisation::delay(BaseType *inFrame, int inLength, BaseType *outFrame, int outLength, int delay)
{
    int length = std::min(inLength, outLength);
    int outstart = std::min(0, outLength - inLength);
    delay = std::min(delay, length);
    delay = std::max(delay, -length);

    wipp::setZeros(outFrame, outLength);

    if (delay < 0)
	wipp::copyBuffer(&inFrame[-delay], &outFrame[outstart],       length + delay);
    else
	wipp::copyBuffer(inFrame,          &outFrame[delay+outstart], length - delay);
}

void TemporalGCCBinauralLocalisation::crossCorrelation(BaseType *leftFrame, int leftLength, BaseType *rightFrame, int rightLength, BaseType *outFrame, int outLength)
{
    wipp::setZeros(outFrame, outLength);
    wipp::cross_corr(leftFrame, leftLength, rightFrame,rightLength, outFrame, outLength, -outLength/2);
    //    ippsCrossCorr_64f(leftFrame, leftLength, rightFrame, rightLength, outFrame, outLength, -outLength/2);
    BaseType stdDevL = 1, stdDevR = 1;
    wipp::stddev(leftFrame, leftLength, &stdDevL);
    wipp::stddev(rightFrame, rightLength, &stdDevR);
    wipp::divC(stdDevL*stdDevR, outFrame, outLength);
}

double TemporalGCCBinauralLocalisation::correlation(BaseType *leftFrame, BaseType *rightFrame, int length)
{
    BaseType aux[length];
    BaseType mean;
    wipp::mult(leftFrame, rightFrame, aux, length);
    wipp::mean(aux, length, &mean);
    return mean;
}

double TemporalGCCBinauralLocalisation::integration(BaseType *inFrame, int length)
{
    BaseType sum;
    wipp::abs(inFrame, length);
    wipp::sum(inFrame, length, &sum);
    return sum;
}

void TemporalGCCBinauralLocalisation::maxminNormalisation(BaseType *vector, int length)
{
    BaseType  max, min;
    wipp::min(vector, length, &min);
    wipp::addC(min, vector, length);
    wipp::max(vector, length, &max);
    wipp::divC(max, vector, length);
}

//
// ((d*SR)/c)* cos(degrees*(PI/180))
// d*SR/c = microphoneDistance*sampleRate/speedOfSound
// d*SR/c = distance between microphones in samples,
// corresponds to the number of delays applied to the signals: _ndelays
//
int TemporalGCCBinauralLocalisation::degrees2Samples(double degrees, int sampleRate)
{
    return _microphoneDistance*cos(degrees*M_PI/180)*sampleRate/getSpeedOfSound();
}

//
// acos(2*(delay- (dist-1))/(dist-1)) * 180/PI
//
// dist is the distance between microphones in samples and is calculated as follows:
// dist = microphoneDistance*samplesRate/speedOfSound = (d*SR/c)
// dist is equivalent to _ndelays (the number of delays applied to the signal)
//
// For _ndelays 8
// Degrees for each delay: 90.0000   45.5847   25.3769    8.2132   -8.2132  -25.3769  -45.5847  -90.0000
// For _ndelays 10
// Degrees for each delay: 90.0000   51.0576   33.7490   19.4712    6.3794   -6.3794  -19.4712  -33.7490  -51.0576  -90.0000
// For _ndelays 17
// Degrees for each delay: 90.0000   58.9973   45.5847   34.8499   25.3769   16.6015    8.2132         0   -8.2132  -16.6015  -25.3769  -34.8499  -45.5847  -58.9973  -90.0000
//
// For a distance of ~6cm between microphoens (0.06m)
// if sampling rate is 44100 ndealys = 8
// if sampling rate is 48000 ndelays = 8
// if sampling rate is 96000 ndelays = 17 <-- desired option
//
double TemporalGCCBinauralLocalisation::samples2Degrees(int delay, int ndelays)
{
    delay = std::min(delay, ndelays-1);
    delay = std::max(delay, -(ndelays-1));
    double shift = static_cast<double>(ndelays-1)/2;
    double degrees = acos(2*(static_cast<double>(delay)-shift)/(ndelays-1))*180/M_PI;
    return degrees;
}

BaseType TemporalGCCBinauralLocalisation::setPowerFloor(std::vector<double*> &analysisFrames, int analysisLength,
							int nchannels, int sampleRate)
{
    int neededSamples = _durationToEstimatePowerFloor*sampleRate;
    analysisLength = std::min(neededSamples - _samplesConsumedForNoise, analysisLength);
    double power = dsp::SignalPower::power(analysisFrames, analysisLength)*analysisLength;
    _powerFloor += power;
    _samplesConsumedForNoise += analysisLength;
    if (_samplesConsumedForNoise >= neededSamples)
    {
	_noiseEstimated = true;
	_powerFloor /= _samplesConsumedForNoise;
	_powerFloor = 0.15*(100 - _powerFloor) + _powerFloor;  // 100 is ~20*log10(32000)
	DEBUG_STREAM("POWER floor: " << _powerFloor);
    }
    return _powerFloor;
}



//----------------- Freq GCC Binaural Localisation ---------------------------------------------------------------

FreqGCCBinauralLocalisation::FreqGCCBinauralLocalisation(int sampleRate, ArrayDescription microphonePositions, bool usePowerFloor) :
    SoundLocalisationImpl(microphonePositions),
    STFTAnalysis(2, calculateOrderFromSampleRate(sampleRate, _frameRate)),
    _corrMemoryFactor(0),
    _doaMemoryFactor(0),
    _microphoneDistance(microphonePositions.distance(0,1)),
    _silenceFramesCounter(0),
    _sampleRate(sampleRate),
    _doaStep(3*M_PI/180),
    _numSteps(round(M_PI/_doaStep) + 1),
    _usePowerFloor(usePowerFloor),
    _gcc(getAnalysisLength()/2, _gcc.ONESIDEDFFT)
{
    if (_microphonePositions.size() != 2)
    {
	WARN_STREAM("The number of microphones in ArrayDescription is different from 2.");
    }
    _currentDOA.reset(new BaseType[1]);
    _prob.reset(new BaseType[1]);
    _currentDOA[0] = 0;
    _prob[0] = -1;

    _sampledDOAs.reset(new uint16_t[_numSteps]);
    _magnitude.reset(new BaseType[getAnalysisLength()/2]);
    _power.reset(new BaseType[getAnalysisLength()/2]);
    _correlations.reset(new BaseTypeC[_numSteps]);
    _correlationsReal.reset(new BaseType[_numSteps]);
    _mixedChannel.reset(new BaseTypeC[getAnalysisLength()]);
    _samplesDelay.reset(new BaseType[_numSteps]);
    _prevCorrelationsReal.reset(new BaseType[_numSteps]);
    _triangle.reset(new BaseType[_numSteps]);

    //    ippsVectorSlope_16u(_sampledDOAs.get(), _numSteps, -M_PI_2, _doaStep);
    wipp::ramp(_sampledDOAs.get(), _numSteps, -M_PI, _doaStep);

    wipp::setZeros(_prevCorrelationsReal.get(), _numSteps);

    double phase = 3*M_PI_2;
    //    ippsTriangle_Direct_64f(_triangle.get(), _numSteps, 0.001, 1/(2*static_cast<double>(_numSteps)), 0, &phase);
    wipp::triangle(_triangle.get(), _numSteps, 2*static_cast<double>(_numSteps), phase);


    // Generate delay grid according to DOA grid
    for (int i=0; i < _numSteps; i++)
    {
	_samplesDelay[i] = doaToDelayFarFieldSamples(doaIdx2angle(i, _doaStep), _microphoneDistance, _sampleRate);
	TRACE_STREAM("Angle grid(" << i << "): " << toDegrees(doaIdx2angle(i, _doaStep)));
    }

    DEBUG_STREAM("MD: " << _microphoneDistance << " DOA step: " << toDegrees(_doaStep));

    _gcc.precomputeTauMatrix(_samplesDelay.get(), _numSteps, getAnalysisLength()/2, _gcc.ONESIDEDFFT);
}

class null_deleter
{
    null_deleter()
    {
      return;
    }
  public:
    ~null_deleter()
    {
      return;
    }
    void operator()(BaseType* v)
    {
      return;
    }
    friend class FreqGCCBinauralLocalisation;
};

BaseType FreqGCCBinauralLocalisation::setPowerFloor(std::vector<BaseType*> &analysisFrames, int analysisLength, int nchannels, int sampleRate)
{
  SignalVector af;
  null_deleter deleter;
  for (size_t i = 0; i < analysisFrames.size(); ++i)
  {
    af.push_back(boost::shared_array<double>(analysisFrames[i], deleter));
  }
  setPowerFloor(af, analysisLength, nchannels, sampleRate);
}

BaseType FreqGCCBinauralLocalisation::setPowerFloor(SignalVector &analysisFrames, int analysisLength, int nchannels, int sampleRate)
{
    int neededSamples = _durationToEstimatePowerFloor*sampleRate;
    double power = dsp::SignalPower::power(analysisFrames, analysisLength)*(analysisLength-2);
    _powerFloor += power + 1e-10;
    _samplesConsumedForNoise += (analysisLength-2); // length is one-sided FFT length here.
    if (_samplesConsumedForNoise >= neededSamples)
    {
	_noiseEstimated = true;
	if (_samplesConsumedForNoise > 0)
	    _powerFloor /= _samplesConsumedForNoise;
	DEBUG_STREAM("POWER floor: " << _powerFloor << " " << _samplesConsumedForNoise);
	_powerFloor = 10*log10(_powerFloor) + _noiseMarginDB;
	DEBUG_STREAM("POWER floor: " << _powerFloor);
    }

    return _powerFloor;
}

void FreqGCCBinauralLocalisation::processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
							 std::vector<double*> &dataChannels, int dataLength)
{

    if (!_ptrCallback)
    {
	ERROR_STREAM("I am not computing binaural localisation because no callback has been set.");
	return;
    }

    // Calculate the phase of the FFT of the window using the GCC function
    // to obtain the DOA estimation
    BaseTypeC *left = reinterpret_cast<BaseTypeC*>(analysisFrames[0]);
    BaseTypeC *right= reinterpret_cast<BaseTypeC*>(analysisFrames[1]);

    BaseType DOA = 0;
    BaseType power = 0;
    BaseType max = 0;

    size_t idx = 0;
    int complexAnalysisLength = analysisLength/2;

    // Getting signal power, and setting as floor if the firsts signal frames are processed.
    if (!_noiseEstimated)
	power = setPowerFloor(analysisFrames, analysisLength, 2, _sampleRate);
    else
	power = dsp::SignalPower::FFTLogPower(analysisFrames, analysisLength);

    if(power > _powerFloor || !_usePowerFloor)
    {

	DEBUG_STREAM( power << " > " << _powerFloor);
	_gcc.calculateCorrelationsForPrecomputedTauMatrix(reinterpret_cast<dsp::Complex*>(left),
							  reinterpret_cast<dsp::Complex*>(right),
							  reinterpret_cast<dsp::Complex*>(_correlations.get()),
							  complexAnalysisLength,
							  _numSteps, _gcc.ONESIDEDFFT);

	wipp::real(reinterpret_cast<wipp::wipp_complex_t*>(_correlations.get()), _correlationsReal.get(), _numSteps);
	wipp::multC(1-_corrMemoryFactor, _correlationsReal.get(), _numSteps);
	wipp::multC(_corrMemoryFactor, _prevCorrelationsReal.get(), _numSteps);
	wipp::add(_prevCorrelationsReal.get(), _correlationsReal.get(), _numSteps);
	wipp::copyBuffer(_correlationsReal.get(), _prevCorrelationsReal.get(), _numSteps);

	//          ippsAdd_64f_I(_triangle.get(), _correlationsReal.get(), _numSteps);



	setProbability(_currentDOA.get(), _prob.get(), 1);

#ifdef USE_PARTICLE_FILTER
	if (!_particleFilter)
	{
	    wipp::maxidx(_correlationsReal.get(), _numSteps, &max, &idx);
	    DOA = doaIdx2angle(idx, _doaStep);
	    _predictionModel.reset(new SoundLocalisationPredicitonModel(DOA, 0));
	    ++_sourceCounter;
	    _particleFilter.reset(new particle_filter_t(DOA, 500, _sourceCounter, std::make_pair(-M_PI,M_PI),
							_observationModel.get(), _predictionModel.get(),
							_resamplingModel.get()));

	    //		_particleFilter.reset(new dsp::ParticleFilter<double, int, >(DOA, 500, _sourceCounter,
	    //									     std::make_pair(-M_PI_2, M_PI_2),
	    //									     _observationModel, _predictionModel));
	    DEBUG_STREAM("Started following source: " << _particleFilter->getId());
	}

	_currentDOA[0] = _particleFilter->updateFilter();


#ifdef PLOT_DEBUG

	static int plotcount = 0;
	if (plotcount % 1 == 0)
	{
	    std::ostringstream oss;
	    oss << "Source " << _sourceCounter;
	    BasicParticleSet<double> particles = _particleFilter->getParticles();
	    BasicParticleSet<double> weights = _particleFilter->getWeights();
	    static Gnuplot plot;
	    plot.reset_all();
	    plot.set_xrange(0,0.2);
	    plot.set_yrange(-0.1,0.1);
	    plot.set_style("points");
	    plot.set_polar();
	    plot.set_title(oss.str());
	    plot.plot_xy(particles.get(), weights.get(), particles.size());
	    BaseType a = 0.1;
	    plot.set_style("impulses");
	    plot.plot_xy(_currentDOA.get(), &a, 1);
	    sys::Timer::sleepMs(50);
	}
	++plotcount;
#endif

#else
	ippsMaxIndx_64f(_correlationsReal.get(), _numSteps, &max, &idx);
	DOA = doaIdx2angle(idx, _doaStep);
	_currentDOA[0] = _doaMemoryFactor * _currentDOA[0] + (1-_doaMemoryFactor)*DOA;
#ifdef PLOT_DEBUG
	static Gnuplot plot;
	plot.reset_all();
	//            plot.set_style("lines");
	//            plot.plot_x(_correlationsReal.get(), _numSteps);
	plot.set_style("impulses");
	plot.set_polar();
	plot.set_xrange(0,0.2);
	plot.set_yrange(-0.1,0.1);
	//            plot.plot_xy(&idx, &max, 1);
	BaseType a = 0.1;
	plot.plot_xy(_currentDOA.get(), &a, 1);
	sys::Timer::sleepMs(50);
#endif
#endif

	_ptrCallback->setDOA(toDegrees(_currentDOA,1), _prob, power,1);

	_corrMemoryFactor = _maxCorrMemoryFactor;
	_doaMemoryFactor = _maxDoaMemoryFactor;
	_silenceFramesCounter = 0;

    }
    else
    {
	if (_noiseEstimated)
	{
	    // Memory factors goes to zero exponentially as silence time increases.
	    int secondsToDecay = 3;
	    int windowsToDecay = secondsToDecay*_sampleRate/(analysisLength/2-1);

	    if (_silenceFramesCounter < windowsToDecay)
	    {
		// Still needs work
		//                _corrMemoryFactor = _maxCorrMemoryFactor * exp(-(static_cast<double>(_silenceFramesCounter++))/(2*windowsToDecay));
		//                _doaMemoryFactor = _maxDoaMemoryFactor * exp(-(static_cast<double>(_silenceFramesCounter++))/(2*windowsToDecay));
		_corrMemoryFactor = _maxCorrMemoryFactor;
		_doaMemoryFactor = _maxDoaMemoryFactor;
		//                DEBUG_STREAM("corr mem factor: " << _corrMemoryFactor << " doa mem factor: " << _doaMemoryFactor << " " << _silenceFramesCounter << " " << windowsToDecay);
		if (_particleFilter)
		{
		    _currentDOA[0] = _particleFilter->updateFilter();
		    _ptrCallback->setDOA(toDegrees(_currentDOA,1), _prob, power,1);
		}

	    }
	    else
	    {
		if (_particleFilter)
		    DEBUG_STREAM("Stopped following source: " << _particleFilter->getId());
		_particleFilter.reset();
		_corrMemoryFactor = 0;
		_doaMemoryFactor = 0;
	    }
	    ++ _silenceFramesCounter;
	}
    }

    TRACE_STREAM("Idx: " << idx << " DOA: " << toDegrees(doaIdx2angle(idx, _numSteps)) << " samples "
		 << doaToDelayFarField(doaIdx2angle(idx, _numSteps), _microphoneDistance)*_sampleRate
		 << ", C: " << max
		 );
}

void FreqGCCBinauralLocalisation::setProbability(const double *doas, double *probs, int size)
{
    double sum;
    double min;
    double angle;
    double prevdoa, nextdoa;

    double probDistr[_numSteps];
    int idx;
    double p;
    double prevcorr, nextcorr;

    // To avoid a warning for no initialisation.
    prevcorr = nextcorr = prevdoa = nextdoa = 0;

    wipp::copyBuffer(_correlationsReal.get(), probDistr, _numSteps);
    wipp::min(probDistr, _numSteps, &min);
    wipp::sum(probDistr, _numSteps, &sum);
    sum -= min*_numSteps; // equivalent to add min to all the values in correlation.


    // here I interpolate beetween the closests values
    // of the sampled space.
    // I do this to avoid having multiple particle in
    // the particle filter with the same weight in the case
    // were they are very close.
    for (int i = 0; i < size; ++i)
    {

	idx = angle2DOAidx(doas[i], _doaStep);
	angle = doaIdx2angle(idx, _doaStep);

	if (0 < idx && idx < (_numSteps -1))
	{
	    if (angle > doas[i] && idx > 0)
	    {
		prevcorr =probDistr[idx-1];
		prevdoa = doaIdx2angle(idx-1, _doaStep);
		nextcorr = probDistr[idx];
		nextdoa = angle;
	    }
	    else if (angle <= doas[i] && idx < (_numSteps -1))
	    {
		prevcorr = probDistr[idx];
		prevdoa = angle;
		nextcorr = probDistr[idx+1];
		nextdoa = doaIdx2angle(idx+1, _doaStep);
	    }

	    double slope = (nextcorr-prevcorr)/(nextdoa-prevdoa);
	    p = slope*(doas[i] - prevdoa) + prevcorr;
	}
	else
	{
	    p = probDistr[idx];
	}

	probs[i] = 0;
	if (sum > 0)
	    probs[i] = (p - min)/sum;
	probs[i] = (probs[i] < 0.01) ? 0 : probs[i];
    }
}



}

