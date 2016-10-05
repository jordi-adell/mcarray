/*
* SteeringBeamforming.cpp
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
#include <mcarray/SteeringBeamforming.h>
#include <mcarray/mcalogger.h>

#include <dspone/algorithm/gralCrossCorrelation.h>
#include <dspone/filter/MedianFilter.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>

namespace mca {

SteeringBeamforming::SteeringBeamforming(int sampleRate, ArrayDescription microphonePositions, int fftCCSLength, unsigned int nchannels) :
  _sampleRate(sampleRate),
  _fftCCSLength(fftCCSLength),
  _complexFFTCCSLength(fftCCSLength/2),
  _nchannels(nchannels),
  _doaStep(5*M_PI/180),
  _numSteps(round(M_PI/_doaStep) + 1),
  _microphonePositions(microphonePositions)
{
  allocate();
  generateLookupTable();
}

void SteeringBeamforming::allocate()
{
  _energyInDOA.reset(new BaseType[_numSteps]);
  _prevEnergyInDOA.reset(new BaseType[_numSteps]);
  wipp::setZeros(_prevEnergyInDOA.get(), _numSteps);
  _firstDerivative.reset(new BaseType[_numSteps-1]);
  _secondDerivative.reset(new BaseType[_numSteps-2]);
  _complexCorrelation.reset(new BaseTypeC[_numSteps]);
}

void SteeringBeamforming::generateLookupTable()
{
  // Calculate the delays (one per DOA) for each pair of microphones (i, j).
  // Precompute delay matrices.

  for (unsigned int i = 0; i < _nchannels; ++i)
  {
    for (unsigned int j = i+1; j < _nchannels; ++j)
    {
      double distance = _microphonePositions.distance(i,j);

      SignalPtr delaysForMicroPair;
      delaysForMicroPair.reset(new BaseType[_numSteps]);
      for (int doa = 0; doa < _numSteps; ++doa)
      {
	delaysForMicroPair[doa] = doaToDelayFarFieldSamples(doaIdx2angle(doa, _doaStep), distance, _sampleRate);

	TRACE_STREAM("Micro Pair [" << i << ", " << j << "]. Delay[doa=" <<
		     toDegrees(doaIdx2angle(doa, _doaStep)) << "]: " << delaysForMicroPair[doa]);
      }

      // One correlation vector for each micro pair. Each position of these vectors will contain the
      // cross-correlation for one delay (DOA).
      _correlations.push_back(SignalPtr(new BaseType[_numSteps]));

      // For each micro pair, create a GCC object and precompute the delay matrix according to 'delaysForMicroPair'
      _gcc.push_back(std::shared_ptr<dsp::GeneralisedCrossCorrelation>
		     (new dsp::GeneralisedCrossCorrelation(_complexFFTCCSLength, dsp::GeneralisedCrossCorrelation::ONESIDEDFFT)));

      _gcc.back()->precomputeTauMatrix(delaysForMicroPair.get(), _numSteps, _complexFFTCCSLength,
				       dsp::GeneralisedCrossCorrelation::ONESIDEDFFT);

      // row 'k' of _microPairIdx contains the indexs of the two microphones of the micro pair 'k'.
      _microPairIdx.push_back({i, j});
    }
  }
}

void SteeringBeamforming::processFrame(const SignalVector &analysisFrames, SignalPtr DOA, SignalPtr prob, int numOfSources, SignalVector &wienerCoefs)
{
  computeCorrelations(analysisFrames, wienerCoefs);

  computeEnergyInDOA();
  selectDOA(DOA, prob, numOfSources);
}

void SteeringBeamforming::computeCorrelations(const SignalVector &analysisFrames, SignalVector &wienerCoefs)
{
  // Compute the correlation for each micro pair.
  for (unsigned int pairIdx = 0; pairIdx < _correlations.size(); ++pairIdx)
  {
    // Take the frames of the micro pair 'pairIdx' using the information in '_microPairIdx'.
    dsp::Complex *frameA = reinterpret_cast<dsp::Complex*>(analysisFrames[_microPairIdx[pairIdx][0]].get());
    dsp::Complex *frameB = reinterpret_cast<dsp::Complex*>(analysisFrames[_microPairIdx[pairIdx][1]].get());

    // Compute the correlations using the precomputed delays matrix.
    // Wiener coefs are not used yet, but are accessible from here.
    _gcc[pairIdx]->calculateCorrelationsForPrecomputedTauMatrix(frameA, frameB,
								reinterpret_cast<dsp::Complex*>(_complexCorrelation.get()),
								_complexFFTCCSLength,
								_numSteps,
								dsp::GeneralisedCrossCorrelation::ONESIDEDFFT);

    // Store the real part of the correlation in _correlations vector.
    wipp::real(reinterpret_cast<wipp::wipp_complex_t*>(_complexCorrelation.get()), _correlations[pairIdx].get(), _numSteps);

    for (int i = 0; i < _numSteps; ++i)
    {
      TRACE_STREAM("Mics Pair [" << _microPairIdx[pairIdx][0] << ", " << _microPairIdx[pairIdx][1] << "]. Corr[doa="
							      << toDegrees(doaIdx2angle(i, _doaStep)) << "]: " << _correlations[pairIdx][i]);
    }
  }
}

void SteeringBeamforming::computeEnergyInDOA()
{
  wipp::multC(_energyMemoryFactor, _prevEnergyInDOA.get(), _energyInDOA.get(), _numSteps);

  // Sum the correlations of all micro pairs (each position corresponds to a DOA).
  for (unsigned int pairIdx = 0; pairIdx < _correlations.size(); ++pairIdx)
  {
    wipp::multC(1-_energyMemoryFactor, _correlations[pairIdx].get(), _numSteps);
    wipp::add(_correlations[pairIdx].get(), _energyInDOA.get(), _numSteps);
  }

  wipp::copyBuffer(_energyInDOA.get(), _prevEnergyInDOA.get(), _numSteps);
}

void SteeringBeamforming::selectDOA(SignalPtr DOA, SignalPtr prob, int numOfSources)
{

  double max;
  size_t maxIdx;
  int numPairs = _correlations.size();
  const double minEnergyInDOA = -15*numPairs; // Min possible value of energyInDOA (ad-hoc).

  // Min-max normalization
  wipp::subC(minEnergyInDOA, _energyInDOA.get(), _numSteps);
  wipp::divC(-2*minEnergyInDOA, _energyInDOA.get(), _numSteps);

  // I derivate the energy estimation  y[n] = x[n] - x[n-1]
  wipp::sub(_energyInDOA.get(), &(_energyInDOA[1]), _firstDerivative.get(), _numSteps-1);
  // I obtain the "sign": 0 for positive values of first derivative, 1 for negative values
  wipp::threshold_lt_gt(_firstDerivative.get(), _numSteps-1, 0.0, 1.0, 0.0, 0.0);
  // Apply median filter to avoid false detections.
  // ippsFilterMedian_64f_I(_firstDerivative.get(), _numSteps-1, 3);
  wipp::median_filter(_firstDerivative.get(), _filteredFirstDerivative.get(), _numSteps-1, 3);

  // I derivate again y[n] = x[n] - x[n-1]. Positive peaks indicates local maximums in energyInDOA.
  // Negative peaks indicates local minimums in energyInDOA.
  wipp::sub(_firstDerivative.get(), &(_firstDerivative[1]), _secondDerivative.get(), _numSteps-2);
  // I multiply the peaks obtained in previous line by _energyInDOA, in order to assign a weight to
  // each one according to the energy.
  wipp::mult(&(_energyInDOA[1]), _secondDerivative.get(), _numSteps-2);

  //          Gnuplot plot("lines");
  //          plot.reset_plot();
  //          plot.plot_x(&(_energyInDOA[1]), _numSteps-2, "energyInDOA");
  //          plot.plot_x(&(_firstDerivative[1]), _numSteps-2, "firstDerivative");
  //          plot.plot_x(_secondDerivative.get(), _numSteps-2, "secondDerivative");
  //          util::Timer::sleepMs(300);

  // Find "numOfSources" maximums in _secondDerivative. Peaks with more energy will be found first.
  // Negative peaks (minimums) will not be found (0 > negative_peak). When all positive peaks are found,
  // then prob = 0.
  for (int s = 0; s < numOfSources; ++s)
  {
    //    ippsMaxIndx_64f(_secondDerivative.get(), _numSteps-2, &max, &maxIdx);
    wipp::maxidx(_secondDerivative.get(), _numSteps-2, &max, &maxIdx);
    _secondDerivative[maxIdx] = 0;

    DOA[s] = doaIdx2angle(maxIdx+1, _doaStep);
    //@TODO: Estimate the probability in a reliable manner.
    prob[s] = max;
  }
}


}

