/*
* SourceSeparationAndLocalisation.h
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

#ifndef __SOURCESEPARATIONANDLOCALISATION_H_
#define __SOURCESEPARATIONANDLOCALISATION_H_


#include <mcarray/BinauralLocalisation.h>
#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/BeamformingSeparationAndLocalistaion.h>


#include <dspone/rt/ShortTimeProcess.h>
#include <dspone/rt/ShortTimeFourierTransform.h>
#include <dspone/algorithm/gralCrossCorrelation.h>

#include <memory>

namespace mca
{


class SourceSeparationAndLocalisation : public dsp::STFT
{

    public:

	SourceSeparationAndLocalisation(int sampleRate, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor=true);

	virtual ~SourceSeparationAndLocalisation(){}

	/**
	   * @brief setCallback
	   * @param callback
	   */
	void setCallback(LocalisationCallback &callback);
	void setCallback(LocalisationCallback *callback);

    private:

	static constexpr float _frameRate = 0.025;  /**< related with the length of the frame to work with, in  seconds  */
	const int _sampleRate; /**< sample rate of the signals to be processed */
	SignalVector _denoisedFrames; /**< frames after noise reduction. */
	std::unique_ptr<BeamformingSeparationAndLocalisation> _impl; /**< pointer to the implementation of beamforming localisation and separation */
	//	NoiseReductionMCRAWiener _noiseReduction; /**< noise reduction object. */
	//	const int _wienerFilterLength; /**< number coeficients of wiener filter. */
	SignalVector _wienerCoefs; /**< vector containing the wiener coeficients. */

	virtual void processParametrisation(SignalVector &analysisFrames, int analysisLength,
					    std::vector<double*> &dataChannels, int dataLength);
	virtual void processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
					    std::vector<double*> &dataChannels, int dataLength);

	/**
	   * @brief allocate Allocates memory.
	   */
	void allocate();
};



//class SourceLocalisation : public dsp::STFTAnalysis
//{

//    public:

//	SourceLocalisation(int sampleRate, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor=true);

//	virtual ~SourceLocalisation(){}

//	/**
//	   * @brief setCallback
//	   * @param callback
//	   */
//	void setCallback(LocalisationCallback &callback);
//	void setCallback(LocalisationCallback *callback);

//    private:

//	static constexpr float _frameRate = 0.025;  /**< related with the length of the frame to work with, in  seconds  */
//	int _sampleRate; /**< sample rate of the signals to be processed */
//	int _wienerFilterLength; /**< number coeficients of wiener filter. */
//	SignalVector _wienerCoefs; /**< vector containing the wiener coeficients. */
//	SignalVector _denoisedFrames; /**< frames after noise reduction. */
//	std::unique_ptr<BeamformingSeparationAndLocalisation> _impl; /**< pointer to the implementation of beamforming localisation and separation */
//	//	NoiseReductionMCRAWiener _noiseReduction; /**< noise reduction object. */

//	virtual void processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
//					    std::vector<double*> &dataChannels, int dataLength);

//	/**
//	   * @brief allocate Allocates memory.
//	   */
//	void allocate();

//};

}

#endif // __SOURCESEPARATIONANDLOCALISATION_H
