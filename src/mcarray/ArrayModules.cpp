/*
* ArrayModules.cpp
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
#include <mcarray/ArrayModules.h>
#include <mcarray/FastBinauralMasking.h>
#include <mcarray/BinauralLocalisation.h>
#include <mcarray/SourceSeparationAndLocalisation.h>
#include <mcarray/SoundLocalisationCallback.h>

namespace mca {

SoundLocalisation::SoundLocalisation(int sampleRate, ArrayDescription microphonePositions, LocalisationCallback *callback)
{
    bool usePowerFloor = true;

    if (microphonePositions.size() == 2)
    {
	FreqGCCBinauralLocalisation *loc = nullptr;
	loc = new FreqGCCBinauralLocalisation(sampleRate, microphonePositions, usePowerFloor);
	if (callback)
	{
	    loc->setCallback(callback);
	}
	_impl.reset(dynamic_cast<dsp::ShortTimeAnalysis*>(loc));
    }
    else
    {
	// Carefull, this 512 should not be here !!! I do not know why it is here.
	BeamformingSeparationAndLocalisation *loc;
	loc = new BeamformingSeparationAndLocalisation(sampleRate, 512, microphonePositions, 1, true);
	if (callback)
	{
	    loc->setCallback(callback);
	}
	_impl.reset(dynamic_cast<dsp::ShortTimeAnalysis*>(loc));
    }

}

SoundLocalisation::~SoundLocalisation()
{

}


BinauralMasking::BinauralMasking(int samplerate,
				 ArrayDescription microphones,
				 float lowFreq,
				 float highFreq,
				 MaskingMethod mmethod,
				 MaskingAlg algorithm)
{
    double microDist = 0;
    if (microphones.size() > 1)
    {
	microDist = microphones.distance(0,1);
    }

    _impl.reset(new FastBinauralMasking(samplerate,
					microDist,
					lowFreq,
					highFreq,
					mmethod,
					algorithm));

}

BinauralMasking::~BinauralMasking()
{

}


}
