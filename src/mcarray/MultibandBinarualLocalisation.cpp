/*
* MultibandBinarualLocalisation.cpp
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
#include <mcarray/MultibandBinarualLocalisation.h>
#include <mcarray/SoundLocalisationCallback.h>

#include <dspone/algorithm/signalPower.h>

#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>
#include <wipp/wipputils.h>

namespace mca{

MultibandBinarualLocalisation::MultibandBinarualLocalisation(int sampleRate, ArrayDescription microphonePositions, int nbins, bool usePowerFloor) :
    SoundLocalisationImpl(microphonePositions),
    dsp::SubBandSTFTAnalysis(nbins,
			     sampleRate,
			     calculateOrderFromSampleRate(sampleRate, _frameRate),
			     2,
			     100,
			     maxFreqForSpatialAliasing(microphonePositions.distance(0,1)),
			     dsp::SubBandSTFT::LINEAR),
    _microphoneDistance(microphonePositions.distance(0,1)),
    _doaStep(5*M_PI/180),  // 5 degree resolution
    _numSteps(floor(M_PI/_doaStep) + 1),
    _sampleRate(sampleRate),
    _numberOfBins(getNumberOfBins()),
    _minFreq(100.0F/_sampleRate),
    _maxFreq(maxFreqForSpatialAliasing(_microphoneDistance)/_sampleRate),
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

    _binDOAs.reset(new BaseType[_numberOfBins]);
    _energyInDOA.reset(new BaseType[_numSteps]);
    _energies.reset(new BaseType[_numberOfBins]);
    _magnitude.reset(new BaseType[getAnalysisLength()/2]);
    _power.reset(new BaseType[getAnalysisLength()/2]);
    _correlations.reset(new BaseTypeC[_numSteps]);
    _correlationsReal.reset(new BaseType[_numSteps]);
    _mixedChannel.reset(new BaseTypeC[getAnalysisLength()]);
    _samplesDelay.reset(new BaseType[_numSteps]);

    _triangle.reset(new BaseType[_numSteps]);
    double phase = 3*M_PI_2;
    //    ippsTriangle_Direct_64f(_triangle.get(), _numSteps, 0.001, 1/(2*static_cast<double>(_numSteps)), 0, &phase);
    wipp::triangle(_triangle.get(), _numSteps, 2*static_cast<double>(_numSteps), phase);
    // Generate delay grid according to DOA grid
    for (int i=0; i < _numSteps; i++) {
	_samplesDelay[i] = doaToDelayFarFieldSamples(doaIdx2angle(i), _microphoneDistance, _sampleRate);
	TRACE_STREAM("Angle grid(" << i << "): " << toDegrees(doaIdx2angle(i)));
    }

    for (int bin = 0; bin < _numberOfBins; ++bin)
    {
	_prevCorrelationsReal.push_back(SignalPtr(new BaseType[_numSteps]));
	wipp::setZeros(_prevCorrelationsReal[bin].get(), _numSteps);
    }

    DEBUG_STREAM("MD: " << _microphoneDistance << " Fmin: " << _filterBank->getBinCenterFrequency(0)*_sampleRate << " Fmax: "
		 << _filterBank->getBinCenterFrequency(_filterBank->getNBins()-1)*_sampleRate
		 << " DOA step: " << toDegrees(_doaStep) << " Nbins: " << _numberOfBins);

    //        _plot.reset(new Gnuplot("lines"));
}


BaseType MultibandBinarualLocalisation::setPowerFloor(std::vector<double*> &analysisFrames, int analysisLength, int nchannels, int sampleRate)
{
    int neededSamples = _durationToEstimatePowerFloor*sampleRate;
    double power = dsp::SignalPower::FFTPower(analysisFrames, analysisLength)*(2*analysisLength-2);
    _powerFloor += power;
    _samplesConsumedForNoise += (2*analysisLength-2); // length is one-sided FFT length here.
    if (_samplesConsumedForNoise >= neededSamples)
    {
	_noiseEstimated = true;
	_powerFloor /= _samplesConsumedForNoise;
	_powerFloor = 10*log10(_powerFloor) + _noiseMarginDB;

	DEBUG_STREAM("POWER floor: " << _powerFloor);
    }

    return _powerFloor;
}

void MultibandBinarualLocalisation::processSetup(std::vector<double *> &analysisFrames, int analysisLength,
						 std::vector<double *> &dataChannels, int dataLength)
{
    wipp::setZeros(_energyInDOA.get(), _numSteps);
    wipp::setZeros(_binDOAs.get(),     _numberOfBins);
    wipp::setZeros(_energies.get(),     _numberOfBins);
}

void MultibandBinarualLocalisation::processOneSubband(const SignalVector &analysisFrame, int length, int bin)
{
    // Calculate the phase of the FFT of the window.
    // Calculate the DOA estimation for the current subband, using the GCC function.
    dsp::Complex *left = reinterpret_cast<dsp::Complex*>(analysisFrame[0].get());
    dsp::Complex *right= reinterpret_cast<dsp::Complex*>(analysisFrame[1].get());

    int complexAnalysisLength = length/2;
    BaseType max;
    size_t idx;

    _gcc.calculateCorrelationsForTauVector(left, right, reinterpret_cast<dsp::Complex*>(_correlations.get()), complexAnalysisLength,
					   _samplesDelay.get(), _numSteps, _gcc.ONESIDEDFFT);


    wipp::real(reinterpret_cast<wipp::wipp_complex_t*>(_correlations.get()), _correlationsReal.get(), _numSteps);
    wipp::multC((1-_corrMemoryFactor), _correlationsReal.get(), _numSteps);
    wipp::multC(_corrMemoryFactor, _prevCorrelationsReal[bin].get(), _numSteps);
    wipp::add(_prevCorrelationsReal[bin].get(), _correlationsReal.get(), _numSteps);
    wipp::copyBuffer(_correlationsReal.get(), _prevCorrelationsReal[bin].get(), _numSteps);
    wipp::maxidx(_correlationsReal.get(), _numSteps, &max, &idx);

    _binDOAs[bin] = doaIdx2angle(idx);

    _energies[bin] = dsp::SignalPower::FFTPower(analysisFrame, length);

    _energyInDOA[idx] += _energies[bin];

    TRACE_STREAM("BIN: " << bin << " " << _filterBank->getBinCenterFrequency(bin)*_sampleRate << "Hz, E: "
		 << _energies[bin] << ", Idx: " << idx << " DOA: " << toDegrees(doaIdx2angle(idx)) << " samples " << doaToDelayFarField(doaIdx2angle(idx), _microphoneDistance)*_sampleRate
		 << ", C: " << max
		 );
}

void MultibandBinarualLocalisation::processSumamry(std::vector<double *> &analysisFrames, int analysisLength,
						   std::vector<double *> &dataChannels, int dataLength)
{
    if (_ptrCallback)
    {
	BaseType DOA = 0;
	BaseType power = 0;
	BaseType maxEnergy;
	size_t idx;

	//  I need to estimate noise Floor based on energies, because _analysisFrames do
	// contain frequencies above maxFreq and thus it containes more energy.
	// Or filter the signal in the parent class.

	// Getting signal power, and setting as floor if the first signal frame is processed.
	if (!_noiseEstimated)
	{
	    power = setPowerFloor(analysisFrames, analysisLength/2, 2, _sampleRate);
	}
	else
	{
	    //              ippsMean_64f(_energies.get(), _numberOfBins, &power);
	    power = dsp::SignalPower::FFTPower(analysisFrames, analysisLength);
	}

	// Only if singal power exceeds the floor power, the DOA is obtained
	if (power > _powerFloor || !_usePowerFloor)
	{
	    wipp::sum(_energyInDOA.get(), _numSteps, &(_prob[0]));
	    wipp::maxidx(_energyInDOA.get(), _numSteps, &maxEnergy, &idx);

	    if (_prob[0] != 0)
	    {
		_prob[0] = _energyInDOA[idx]/_prob[0];
	    }

	    //              debugPlot(idx);

	    DOA = doaIdx2angle(idx);

	    _currentDOA[0] = _doaMemoryFactor * _currentDOA[0] + (1-_doaMemoryFactor)*DOA;

	    // This is to change fast to a doa, if we are at 0, wich means default DOA.
	    // However, this might cause problems when someone speaking at front.
	    //              if (-_degreeStep < _currentDOA && _currentDOA < _degreeStep)
	    //                _currentDOA = DOA;
	    //              else
	    //                _currentDOA = _doaMemoryFactor * _currentDOA + (1-_doaMemoryFactor)*DOA;

	    _ptrCallback->setDOA(toDegrees(_currentDOA,1), _prob, power,1);

	}
	else
	{
	    // if the signal power is not high enough DOA is cosidered 0º (strict front)
	    _currentDOA[0] = _currentDOA[0]*_doaMemoryFactorSilence + (1 - _doaMemoryFactorSilence)*0;
	    _prob[0] = -100000;
	    //             _ptrCallback->setDOA(DOA, prob, power);
	}

    }
}

//      void MultibandBinarualLocalisation::debugPlot(int idx)
//      {
//        BaseType max;
//        ippsLn_64f_I(_energies.get(), _numberOfBins);
//        ippsMax_64f(_energies.get(), _numberOfBins, &max);
//        ippsDivC_64f_I(max, _energies.get(), _numberOfBins);

//        //          ippsDivC_64f_I(M_PI, _binDOAs.get(), _numberOfBins);

//        ippsMax_64f(_energyInDOA.get(), _numSteps, &max);
//        ippsDivC_64f_I(max, _energyInDOA.get(), _numSteps);

//        _plot->reset_plot();
//        _plot->set_style("lines");
//        _plot->set_grid();
//        _plot->set_smooth("csplines");
//        _plot->plot_x(_binDOAs.get(), _numberOfBins, "DOA");
//        _plot->plot_x(_energies.get(), _numberOfBins, "E");
//        _plot->plot_x(_energyInDOA.get(), _numSteps, "E/DOA");
//        _plot->plot_x(_triangle.get(), _numSteps, "DOA offset");
//        TRACE_STREAM("nbins: " << _numberOfBins);
//        TRACE_STREAM("steps: " << _numSteps);

//        _plot->set_style("points");
//        int one=1;
//        _plot->unset_smooth();
//        _plot->set_pointsize(5);
//        _plot->plot_xy(&idx, &one , 1, "Selected DOA");
//        idx = angle2DOAidx(_currentDOA[0]);
//        _plot->plot_xy(&idx, &one , 1, "Published DOA");

//        util::Timer::sleepMs(600);
//      }

}

