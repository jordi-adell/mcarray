/*
* BeamformingSeparationAndLocalisation.cpp
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
#include <mcarray/BeamformingSeparationAndLocalistaion.h>
#include <dspone/algorithm/signalPower.h>

#include <wipp/wipputils.h>

namespace mca {

BeamformingSeparationAndLocalisation::BeamformingSeparationAndLocalisation(int sampleRate, int fftCCSLength, ArrayDescription microphonePositions,
									   unsigned int numOfSources, bool usePowerFloor) :
  SoundLocalisationImpl(microphonePositions),
  _nchannels(microphonePositions.size()),
  _sampleRate(sampleRate),
  _fftCCSLength(fftCCSLength),
  _usePowerFloor(usePowerFloor),
  _numOfSources(numOfSources),
  _steeringBeamforming(sampleRate, microphonePositions, fftCCSLength, _nchannels),
  _beamformer(sampleRate, microphonePositions, fftCCSLength, _nchannels)
{
  allocate();
}

void BeamformingSeparationAndLocalisation::allocate()
{
  for (unsigned int c = 0; c < _nchannels; ++c)
    _inputFrames.push_back(SignalPtr(new BaseType[_fftCCSLength]));

  _currentDOA.reset(new BaseType[_numOfSources]);
  _prob.reset(new BaseType[_numOfSources]);

  wipp::setZeros(_currentDOA.get(), _numOfSources);
  wipp::set(-1.0, _prob.get(), _numOfSources);
}

BaseType BeamformingSeparationAndLocalisation::setPowerFloor(SignalVector &analysisFrames, int fftCCSLength, int nchannels, int sampleRate)
{
  int neededSamples = _durationToEstimatePowerFloor*sampleRate;
  double power = dsp::SignalPower::FFTPower(analysisFrames, fftCCSLength)*(fftCCSLength-2);
//      dsp::calculateLinearPowerFFT(analysisFrames, fftCCSLength, nchannels)*(fftCCSLength-2);
  _powerFloor += power;
  _samplesConsumedForNoise += (fftCCSLength-2);
  if (_samplesConsumedForNoise >= neededSamples)
  {
    _noiseEstimated = true;
    _powerFloor /= _samplesConsumedForNoise;
    _powerFloor = 10*log10(_powerFloor) + _noiseMarginDB;

    DEBUG_STREAM("POWER floor: " << _powerFloor);
  }

  return _powerFloor;
}

void BeamformingSeparationAndLocalisation::processFrameLocalisation(SignalVector &analysisFrames, SignalVector &wienerCoefs)
{
  BaseType power;

  // _powerFloor is estimated during an initial number of frames (only if _usePowerFloor==true).
  // Once powerFloor is estimated (_noiseEstimated==true), the power of the current frame is obtained.
  if (!_noiseEstimated && _usePowerFloor)
    power = setPowerFloor(analysisFrames, _fftCCSLength, _nchannels, _sampleRate);
  else
    power = dsp::SignalPower::FFTLogPower(analysisFrames, _fftCCSLength);
  //	dsp::calculateLogPowerFFT(analysisFrames, _fftCCSLength, _nchannels);

  // If the power of the current frame is greater than _powerFloor, then the frame is processed.
  if ((power > _powerFloor) || !_usePowerFloor)
  {
    _steeringBeamforming.processFrame(analysisFrames, _currentDOA, _prob, _numOfSources, wienerCoefs);
    // We publish only the DOA of the first source.
    _ptrCallback->setDOA(toDegrees(_currentDOA, _numOfSources), _prob, power, _numOfSources);
  }
}

void BeamformingSeparationAndLocalisation::processFrameSeparation(SignalVector &inputFrames, SignalVector &outputFrames)
{
  unsigned int c;

  // To prevent incorrect behaivor when inputFrames and outputFrames are the same pointer
  // (inputFrames cannot be modified until all outputFrames are computed)
  for (c = 0; c < _nchannels; ++c)
    wipp::copyBuffer(inputFrames[c].get(), _inputFrames[c].get(), _fftCCSLength);

  // Compute "min(_nchannels, _numOfSources)" outputs and store them in outputFrames.
  for (c = 0; c < std::min(_nchannels, _numOfSources); ++c)
    _beamformer.processFrame(_inputFrames, outputFrames[c], _currentDOA[c]);

  // Set to zero the remaining channels (when _numOfSources < _nchannels).
  for (; c < _nchannels; ++c)
    wipp::setZeros(outputFrames[c].get(), _fftCCSLength);
}



}

