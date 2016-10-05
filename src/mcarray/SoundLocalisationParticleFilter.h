/*
* SoundLocalisationParticleFilter.h
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

#ifndef __SOUNDLOCALISATIONPARTICLEFILTER_H_
#define __SOUNDLOCALISATIONPARTICLEFILTER_H_

#include <mcarray//SoundLocalisationImpl.h>
#include <dspone/pf/ParticleFilter.hpp>
#include <dspone/pf/ParticleSet.hpp>
#include <dspone/pf/PredictionModel.hpp>

namespace mca {

/**
	 * @brief The SoundLocalisationObservationModel class
	 * This class will be used in the particle filter to
	 * get the weights of a set of particles given the the localisation
	 * measurments.
	 */
class SoundLocalisationObservationModel : public dsp::IObservationModel<dsp::ParticleSet<BaseType> >
{

    public:
	/**
	   * @brief SoundLocalisationObservationModel
	   * Constructs and observation model from a SoundLocalisationImpl algorithm.
	   * the function setProbability is used to set the weights of the DOAs.
	   * DOAs a treated as particles in the particle filter.
	   * @param soundLocalisation   a sound localisation algorithm
	   */
	SoundLocalisationObservationModel(SoundLocalisationImpl *soundLocalisation);
	virtual ~SoundLocalisationObservationModel();

	/**
	   * @brief getWeights   returns a weight value for each particle in the particle set
	   * it used the Localisation algorithm to set the probability of each particle.
	   * @param particle
	   * @return
	   */
	virtual dsp::ParticleSet<BaseType> getWeights(const dsp::ParticleSet<BaseType> &particle) const;

	/**
	   * @brief updateModel
	   * This functino is called from the particle filter, to give the oportunity
	   * to get new observations. However, in the design of the sound localisation
	   * the sound localisation algorithm drives the application, so this function is not
	   * doing anything.
	   */
	virtual void updateModel();

    private:
	SoundLocalisationImpl *_soundLocalisation; /**< A pointer of the sound localisation algoritmh **/
};


class SoundLocalisationPredicitonModel : public dsp::PredictionModel<BaseType>
{
    public:
	SoundLocalisationPredicitonModel(BaseType initialState, BaseType initialVelocity);
	virtual ~SoundLocalisationPredicitonModel();
	virtual void update(dsp::ParticleSet<BaseType> &particles);
};


}

#endif // SOUNDLOCALISATIONPARTICLEFILTER_H
