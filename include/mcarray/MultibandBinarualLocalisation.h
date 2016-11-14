/*
* MultibandBinarualLocalisation.h
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

#ifndef __MULTIBANDBINARUALLOCALISATION_H_
#define __MULTIBANDBINARUALLOCALISATION_H_

#include "BinauralLocalisation.h"

#include <mcarray/microhponeArrayHelpers.h>

#include <dspone/filter/FilterBank.h>
#include <dspone/rt/ShortTimeFourierSubBand.h>
#include <dspone/algorithm/gralCrossCorrelation.h>

namespace mca {

class LocalisationCallback;

class MultibandBinarualLocalisation : public SoundLocalisationImpl, public dsp::SubBandSTFTAnalysis
{
    public:
	MultibandBinarualLocalisation(int sampleRate, ArrayDescription microphonePositions, int nbins=15, bool userPowerFloor=1);

    private:
	static constexpr float _frameRate = 0.025;
	static constexpr float _doaMemoryFactor = 0;
	static constexpr float _doaMemoryFactorSilence = 1;
	static constexpr float _corrMemoryFactor = 0.4;
	static constexpr float _noiseMarginDB = 3.0; /**< margin over the noise level to decide that signal is present  */
	const double _microphoneDistance;
	const float _doaStep;
	const int _sampleRate;
	const int _numberOfBins;
	int _numSteps;
	const float _minFreq;
	const float _maxFreq;
	bool _usePowerFloor;
	//	SoundSampleTypeVector _filteredSignals;
	SignalPtr _binDOAs;
	SignalPtr _energyInDOA;
	SignalPtr _energies;
	SignalPtr _magnitude;
	SignalPtr _power;
	SignalCPtr _correlations;
	SignalPtr _correlationsReal;
	SignalPtr _triangle;
	SignalCPtr _mixedChannel;
	SignalPtr _samplesDelay;
	SignalVector _prevCorrelationsReal;

	dsp::GeneralisedCrossCorrelation _gcc;

	virtual void processSetup(std::vector<double *> &analysisFrames, int analysisLength,
				  std::vector<double *> &dataChannels, int dataLength);
	void processOneSubband(const SignalVector &analysisFrame, int length, int bin);
	virtual void processOneSubband(std::vector<double*> &analysisFrame, int length, int bin);

	virtual void processSumamry(std::vector<double *> &analysisFrames, int analysisLength,
				    std::vector<double *> &dataChannels, int dataLength);

	BaseType calculateLinearPower(const BaseTypeC *left, const BaseTypeC *right, int length);
	BaseType calculateLogPower(const BaseTypeC *left, const BaseTypeC *right, int length);
	BaseType setPowerFloor(std::vector<double*> &analysisFrames, int analysisLength, int nchannels, int sampleRate);
	BaseType setPowerFloor(SignalVector &analysisFrames, int analysisLength, int nchannels, int sampleRate);


	inline int angle2DOAidx(float angle) const
	{
	    return (static_cast<int>((angle + M_PI_2)/_doaStep));
	}

	inline float doaIdx2angle(int idx) const
	{
	    float angle = (static_cast<float>(idx)*_doaStep)-M_PI_2;
	    return angle;
	}

	//        void debugPlot(int idx);

};

}



#endif // __MULTIBANDBINARUALLOCALISATION_H
