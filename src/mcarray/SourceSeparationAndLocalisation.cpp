/*
* SourceSeparationAndLocalisation.cpp
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
#include <mcarray/SourceSeparationAndLocalisation.h>
#include <mcarray/SoundLocalisationCallback.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>

#include <math.h>

namespace mca {

SourceSeparationAndLocalisation::SourceSeparationAndLocalisation(int sampleRate, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor) :
  STFT(microphonePositions.size(), calculateOrderFromSampleRate(sampleRate, _frameRate)),
  _sampleRate(sampleRate)
//  _noiseReduction(sampleRate, microphonePositions.size(), _analysisLength),
// _wienerFilterLength(_noiseReduction.getWienerFilterLength())
{
  _impl.reset(new BeamformingSeparationAndLocalisation(sampleRate, getAnalysisLength(), microphonePositions, numOfSources, usePowerFloor));

  for (unsigned int c = 0; c < getNumberOfChannels(); c++)
  {
    _wienerCoefs.push_back(SignalPtr(new BaseType[getAnalysisLength()]));
    _denoisedFrames.push_back(SignalPtr(new BaseType[getAnalysisLength()]));
  }
}

void SourceSeparationAndLocalisation::processParametrisation(
    SignalVector &analysisFrames, int analysisLength,
    std::vector<double*> &dataChannels, int dataLength)
{
  // Apply noise reduction
  //  _noiseReduction.processFrame(_analysisFrames, _denoisedFrames);
  //  _noiseReduction.getWienerCoefs(_wienerCoefs);

  // Sound Localisation
  _impl->processFrameLocalisation(analysisFrames, _wienerCoefs);
  //          _impl->processFrameLocalisation(_denoisedFrames, _wienerCoefs);

  // Source Separation
  _impl->processFrameSeparation(analysisFrames, analysisFrames);
  //          _impl->processFrameSeparation(_denoisedFrames, _analysisFrames);
}

void SourceSeparationAndLocalisation::setCallback(LocalisationCallback *callback)
{
  _impl->setCallback(callback);
}

void SourceSeparationAndLocalisation::setCallback(LocalisationCallback &callback)
{
  _impl->setCallback(callback);
}

//----------------- SourceLocalisation ---------------------------------------------------

//SourceLocalisation::SourceLocalisation(int sampleRate, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor) :
//  STFTAnalysis(microphonePositions.size(), calculateOrderFromSampleRate(sampleRate, _frameRate)),
//  _sampleRate(sampleRate)
////  _noiseReduction(sampleRate, microphonePositions.size(), _analysisLength)
//{
//  _impl.reset(new BeamformingSeparationAndLocalisation(sampleRate, getAnalysisLength(), microphonePositions, numOfSources, usePowerFloor));

//  //_wienerFilterLength = _noiseReduction.getWienerFilterLength();

//  for (unsigned int c = 0; c < getNumberOfChannels(); c++)
//  {
//    _wienerCoefs.push_back(SignalPtr(new BaseType[_wienerFilterLength]));
//    _denoisedFrames.push_back(SignalPtr(new BaseType[getAnalysisLength()]));
//  }
//}


//class fake_deleter
//{
//  public:
//    fake_deleter ()
//    {
//      return;
//    }

//    ~fake_deleter()
//    {
//      return;
//    }

//    void operator ()(double *c)
//    {
//      return;
//    }
//    static fake_deleter del;
//};

////typedef boost::shared_array<double, fake_deleter> non_delete_shared_ptr;


//void SourceLocalisation::processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
//						std::vector<double*> &dataChannels, int dataLength)
//{
//  // Apply noise reduction
//  //  _noiseReduction.processFrame(_analysisFrames, _denoisedFrames);
//  // _noiseReduction.getWienerCoefs(_wienerCoefs);

//  // Sound Localisation

//  SignalVector af;
//  SignalVector dc;

//  for (size_t i = 0; i < analysisFrames.size(); ++i)
//  {
//    af.push_back(boost::shared_array<double>(analysisFrames[i], fake_deleter::del));
//  }

//  for (size_t i = 0; i < dataChannels.size(); ++i)
//  {
//    dc.push_back(boost::shared_array<double>(dataChannels[i], fake_deleter::del));
//  }

//  _impl->processFrameLocalisation(af, _wienerCoefs);
//  //          _impl->processFrameLocalisation(_denoisedFrames, _wienerCoefs);
//}

//void SourceLocalisation::setCallback(LocalisationCallback *callback)
//{
//  _impl->setCallback(callback);
//}


}

