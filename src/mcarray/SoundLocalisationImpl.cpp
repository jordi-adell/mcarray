/*
* SoundLocalisationImpl.cpp
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
#include <mcarray/SoundLocalisationImpl.h>
#include <mcarray/SoundLocalisationCallback.h>
#include "SoundLocalisationParticleFilter.h"

#include <wipp/wipputils.h>

namespace mca {


SoundLocalisationImpl::SoundLocalisationImpl(ArrayDescription microphonePositions) :
  _ptrCallback(NULL),
  _microphonePositions(microphonePositions),
  _powerFloor(0),
  _noiseEstimated(false),
  _samplesConsumedForNoise(0)
{
  _observationModel.reset(new SoundLocalisationObservationModel(this));
  _predictionModel.reset(new SoundLocalisationPredicitonModel(0, 0));
  _resamplingModel.reset(dsp::make_resampling_model<double,double>());
  _sourceCounter=0;
}


void SoundLocalisationImpl::setCallback(LocalisationCallback *callback)
{
  _ptrCallback = callback;
}

void SoundLocalisationImpl::setCallback(LocalisationCallback &callback)
{
  _ptrCallback = &callback;
}

void SoundLocalisationImpl::setProbability(const double* , double *probs, int size)
{
  wipp::set(static_cast<double>(1)/size, probs, size);
}


void SoundLocalisationImpl::setParticleDOAs(std::vector<double> &doas, std::vector<double> &weights) const
{
  if (_particleFilter)
  {
    dsp::ParticleSet<double> pfParticles = _particleFilter->getParticles();
    dsp::ParticleSet<double> pfWeights = _particleFilter->getWeights();

    for (size_t i = 0; i < pfParticles.size() && i < pfWeights.size(); ++i)
    {
      doas.push_back(pfParticles.at(i));
      weights.push_back(pfWeights.at(i));
    }
  }
}


}
