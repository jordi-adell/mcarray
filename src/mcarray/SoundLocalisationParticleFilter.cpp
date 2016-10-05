/*
* SoundLocalisationParticleFilter.cpp
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
#include "SoundLocalisationParticleFilter.h"

#include <dspone/pf/ParticleSet.hpp>
#include <dspone/pf/ParticleFilter.hpp>
#include <dspone/pf/PredictionModel.hpp>

#include <wipp/wippsignal.h>

namespace mca{


// -------- SoundLocalisationObservationModel ------------------------------


SoundLocalisationObservationModel::SoundLocalisationObservationModel(SoundLocalisationImpl *soundLocalisation)  :
    _soundLocalisation(soundLocalisation)
{

}

SoundLocalisationObservationModel::~SoundLocalisationObservationModel()
{

}

dsp::ParticleSet<BaseType> SoundLocalisationObservationModel::getWeights(const dsp::ParticleSet<BaseType> &particle) const
{
    int size = particle.size();
    dsp::ParticleSet<BaseType> weights(size);
    _soundLocalisation->setProbability(particle.get(), weights.get(), size);
    return weights;
}

void SoundLocalisationObservationModel::updateModel()
{

}

// -------- SoundLocalisationPredictionModel -------------------------------

SoundLocalisationPredicitonModel::SoundLocalisationPredicitonModel(BaseType initialState, BaseType initialVelocity) :
    dsp::PredictionModel<BaseType>(initialState, initialVelocity)
{

}

SoundLocalisationPredicitonModel::~SoundLocalisationPredicitonModel()
{

}

void SoundLocalisationPredicitonModel::update(dsp::ParticleSet<BaseType> &particles)
{
    dsp::PredictionModel<BaseType>::update(particles);
    wipp::threshold_lt_gt(particles.get(), particles.size(), -M_PI_2, -M_PI_2, M_PI_2, M_PI_2);
}

}

