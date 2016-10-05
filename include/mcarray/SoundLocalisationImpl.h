/*
* SoundLocalisationImpl.h
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
#ifndef __SOUNDLOCALISATION_H_
#define __SOUNDLOCALISATION_H_

#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/ArrayDescription.h>
#include <mcarray/SoundLocalisationCallback.h>
#include <mcarray/mcadefs.h>

#include <dspone/pf/ParticleFilter.hpp>
#include <dspone/pf/ParticleSet.hpp>
#include <dspone/pf/ResamplingModel.hpp>

namespace mca {

class SoundLocalisationObservationModel;
class SoundLocalisationPredicitonModel;
class LocalisationCallback;

class SoundLocalisationImpl
{
    public:
	/**
	   * @brief SoundLocalisationImpl  constructs the class
	   * @param microphonePositions  Description of the microphone array. Contains the
	   * position of each microphone, in m (resolution of mm is recommended).
	   */
	SoundLocalisationImpl(ArrayDescription microphonePositions);
	virtual ~SoundLocalisationImpl(){}

	/**
	   * @brief setCallback  set the callback object whose setDOA function will be called.
	   * @param callback  a reference to the callback object funtion.
	   */
	void setCallback(LocalisationCallback &callback);
	/**
	   * @brief setCallback  set the callback object whose setDOA function will be called.
	   * @param callback  a reference to the callback object funtion.
	   */
	void setCallback(LocalisationCallback *callback);

	/**
	   * @brief setProbability    set the probability of the provided DOAs
	   * given the current observations. This function will be used
	   * by the observation model to estimate the probabilities of a set of DOAs.
	   * You MUST overload it for each localisation method.
	   * @param doas   DOAs to be estimated its probability
	   * @param probs   probability for each DOA
	   * @param size  length of the vectors.
	   */
	virtual void setProbability(const double *doas, double *probs, int size);

	void setParticleDOAs(std::vector<double> &doas, std::vector<double> &weights) const;

    protected:

	static constexpr double _durationToEstimatePowerFloor = 3; /**< amount of time in seconds used to estimat the noise floor */

	LocalisationCallback* _ptrCallback; /**< pointer to the callback whose setDAO function will be called */
	const ArrayDescription _microphonePositions; /**< array containing the postition of each microphone in the array, in meters */

	SignalPtr _currentDOA; /**< last calculated DOAs, in degrees */
	SignalPtr _prob; /**< Probability associated to each DOA */
	double _powerFloor;  /**< minimum frame power for the setDOA function to be called */
	bool _noiseEstimated; /**< keeps track wether the amount of sampels for noise floor estimation have been consumed */
	int _samplesConsumedForNoise; /**< number of smaples consumend to calculated the noise floor */
	long int _sourceCounter;

	typedef dsp::ParticleFilter<double, int, dsp::ParticleSet<double>, dsp::ParticleSet<double> > particle_filter_t;
	std::unique_ptr<particle_filter_t > _particleFilter;
	std::shared_ptr<SoundLocalisationObservationModel> _observationModel;
	std::shared_ptr<SoundLocalisationPredicitonModel> _predictionModel;
	std::shared_ptr<dsp::ResamplingModel<double, double> > _resamplingModel;

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
	virtual BaseType setPowerFloor(SignalVector &analysisFrames, int analysisLength, int nchannels, int sampleRate) = 0;

};


}


#endif // __SOUNDLOCALISATION_H
