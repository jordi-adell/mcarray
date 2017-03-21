/*
* Beamformer.cpp
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
#include <mcarray/Beamformer.h>
#include <mcarray/microhponeArrayHelpers.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
//#include <wipp/wippstats.h>

//#include <math.h>

namespace mca {

Beamformer::Beamformer(int sampleRate, ArrayDescription microphonePositions, int fftCCSLength, unsigned int nchannels) :
  _sampleRate(sampleRate),
  _fftCCSLength(fftCCSLength),
  _nchannels(nchannels),
  _microphonePositions(microphonePositions)
{
  allocate();
}

void Beamformer::allocate()
{
  _channelSignal.reset(new BaseTypeC[_fftCCSLength/2]);
  _phase.reset(new BaseType[_fftCCSLength/2]);
  _complexRamp.reset(new BaseTypeC[_fftCCSLength/2]);
  _ones.reset(new BaseType[_fftCCSLength/2]);
  wipp::set(1.0, _ones.get(), _fftCCSLength/2);
}

void Beamformer::processFrame(SignalVector &analysisFrames, SignalPtr outputFrame, double DOA)
{
    wipp::setZeros(outputFrame.get(), _fftCCSLength);

    // Sum all channels multiplied by a phase-ramp according to the delay to be applied over each one.
    for (int c = 0; c < _nchannels; ++c)
    {
      //	ippsVectorRamp_64f(_phase.get(),  _fftCCSLength/2, 0, 2*M_PI*_sampleRate/(_fftCCSLength-2)/getSpeedOfSound()*_microphonePositions[c]*cos(DOA+M_PI/2));
	wipp::ramp(_phase.get(), _fftCCSLength/2, 0, 2*M_PI*_sampleRate/(_fftCCSLength-2)/getSpeedOfSound()*_microphonePositions.getX(c)*cos(DOA+M_PI/2));
	wipp::polar2cart(_ones.get(), _phase.get(), reinterpret_cast<wipp::wipp_complex_t*>(_complexRamp.get()), _fftCCSLength/2);
	wipp::mult(reinterpret_cast<wipp::wipp_complex_t*>(analysisFrames[c].get()),
		   reinterpret_cast<wipp::wipp_complex_t*>(_complexRamp.get()),
		   reinterpret_cast<wipp::wipp_complex_t*>(_channelSignal.get()),
		   _fftCCSLength/2);
	wipp::add(reinterpret_cast<wipp::wipp_complex_t*>(_channelSignal.get()),
		  reinterpret_cast<wipp::wipp_complex_t*>(outputFrame.get()),
		  _fftCCSLength/2);
    }

    wipp::divC(_nchannels, outputFrame.get(), _fftCCSLength);
}



}

