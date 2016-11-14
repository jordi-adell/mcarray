/*
* SourceLocalisation.cpp
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
#include <mcarray/SourceLocalisation.h>
#include <mcarray/SoundLocalisationCallback.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>


namespace mca {


class null_deleter
{
    null_deleter()
    {
      return;
    }
  public:
    ~null_deleter()
    {
      return;
    }
    void operator()(BaseType* v)
    {
      return;
    }
    friend class SourceLocalisation;
};

SourceLocalisation::SourceLocalisation(int sampleRate, ArrayDescription microphonePositions, unsigned int numOfSources, bool usePowerFloor) :
  STFTAnalysis(microphonePositions.size(), dsp::ShortTimeProcess::calculateOrderFromSampleRate(sampleRate, _frameRate)),
  _sampleRate(sampleRate)
{
  _impl.reset(new BeamformingSeparationAndLocalisation(sampleRate, getAnalysisLength(), microphonePositions, numOfSources, usePowerFloor));

  for (unsigned int c = 0; c < getNumberOfChannels(); c++)
  {
    _denoisedFrames.push_back(SignalPtr(new BaseType[getAnalysisLength()]));
    _subBandWeights.push_back(SignalPtr(new BaseType[getAnalysisLength()]));
  }
}

void SourceLocalisation::processParametrisation(std::vector<double *> &analysisFrames, int analysisLength,
						std::vector<double *> &dataChannels, int dataLength)
{
  // Apply noise reduction
    //  _noiseReduction.processFrame(_analysisFrames, _denoisedFrames);
    //  _noiseReduction.getWienerCoefs(_wienerCoefs);

  // Sound Localisation
  std::vector<boost::shared_array<double> > af;
  null_deleter deleter;
  for (auto it = analysisFrames.begin(); it != analysisFrames.end(); ++it)
  {
    af.push_back(boost::shared_array<double>(*it, deleter));
  }
  _impl->processFrameLocalisation(af, _subBandWeights);
  //          _impl->processFrameLocalisation(_denoisedFrames, _wienerCoefs);

}

void SourceLocalisation::setCallback(LocalisationCallback *callback)
{
  _impl->setCallback(callback);
}


}

