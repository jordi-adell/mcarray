/*
* BeamformingSeparationAndLocalistaion.h
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
#ifndef __MCA_BEAMFORMING_SEPARATIONANDLOCALISATION_H_
#define __MCA_BEAMFORMING_SEPARATIONANDLOCALISATION_H_


#include <mcarray/SoundLocalisationImpl.h>
#include <mcarray/Beamformer.h>
#include <mcarray/SteeringBeamforming.h>
#include <memory>
#include <vector>

namespace mca
{

class BeamformingSeparationAndLocalisation : public SoundLocalisationImpl
{

    public:

	BeamformingSeparationAndLocalisation(int sampleRate, int fftCCSLength, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor);

	virtual ~BeamformingSeparationAndLocalisation(){}
	void processFrameLocalisation(SignalVector &analysisFrames, SignalVector &wienerCoefs);
	void processFrameSeparation(SignalVector &analysisFrames, SignalVector &outputFrames);

    private:

	const unsigned int _nchannels; /**< number of channels (microphones in array) */
	int _sampleRate; /**< sample rate of the signals to be processed */
	int _fftCCSLength; /**< real length of the buffers containing the FFT in CCS format */
	bool _usePowerFloor; /**< flag used to indicate if power floor is used */
	static constexpr double _noiseMarginDB = 3; /**< noise margin respect to the power floor (in dB). */
	unsigned int _numOfSources; /**< num of sources to search. */
	SignalVector _inputFrames; /**< used to store a copy of inputFrames in processFrameSeparation function. */
	SteeringBeamforming _steeringBeamforming; /**< steering beamforming object.*/
	Beamformer _beamformer; /** beamformer object. */


	/**
	   * @brief setPowerFloor  Sets the power floor and when the necessary amount of samples have been
	   * consumed set the variable __noiseEstimated to true and the variable _powerFloor.
	   * Used calculatePower to obtain the power value.
	   * to the estimated value.
	   * @param analysisFrames  input frames
	   * @param length  length of the input vectors
	   * @param nchannels  number of channels in analysisFrames
	   * @param sampleRate  sampleRate
	   * @return   returns the power calculated.
	   */
	BaseType setPowerFloor(SignalVector &analysisFrames, int fftCCSLength, int nchannels, int sampleRate);

	/**
	   * @brief allocate Allocates memory.
	   */
	void allocate();
};

}


#endif
