/*
* SteeringBeamforming.h
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
#ifndef __MCA_STERING_BEAMFORMING_H_
#define __MCA_STERING_BEAMFORMING_H_

#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/ArrayDescription.h>

#include <memory>

namespace dsp
{
class GeneralisedCrossCorrelation;
}

namespace mca
{


class SteeringBeamforming
{
    public:

	SteeringBeamforming(int sampleRate, ArrayDescription microphonePositions, int fftCCSLength, unsigned int nchannels);
	virtual ~SteeringBeamforming(){}

	/**
	   * @brief processFrame  Processes one frame and computes the DOA estimation and calls the callback with the obtined values.
	   * @param analysisFrames Frames to be processed (FFT in CSS format).
	   * @param DOA  Vector with the DOA estimations, one for each source.
	   * @param prob  Probabiltiy assigned to each DOA.
	   * @param numOfSources  Number of sources to search.
	   * @param wienerCoefs  Wiener coeficients to be used in the GCC (not implemented yet).
	   */
	void processFrame(const SignalVector &analysisFrames, SignalPtr DOA, SignalPtr prob, int numOfSources, SignalVector &wienerCoefs);

    private:

	int _sampleRate; /**< sample rate of the signals to be processed */
	int _fftCCSLength; /**< real length of the buffers containing the FFT in CCS format */
	int _complexFFTCCSLength; /**< complex length of the buffers containing the FFT in CCS format */
	const unsigned int _nchannels; /**< number of channels in input signal (number of microphones) */
	const float _doaStep; /**< step in degrees between DOA values used in the DOA grid. It determines the resolution. */
	const int _numSteps; /**< number of steps in DOA grid (calculated from doaStep). */
	SignalVector32 _wienerCoefs; /**< vector containing the wiener coeficients. */
	ArrayDescription _microphonePositions; /**< description of the array (position of each microphone). */
	SignalCPtr _complexCorrelation; /**< buffer used to compute the correlations. */
	SignalVector _correlations; /**< vector used to store the correlation results for each micro pair. */
	SignalPtr _energyInDOA; /**<  vector that contains the energy in each DOA. */
	SignalPtr _prevEnergyInDOA; /** < used to average with previous energy vector */
	static constexpr float _energyMemoryFactor = 0.8; /**< memory factor to smooth the energy (weigth assigned to the previous energy).  */
	SignalPtr _firstDerivative; /**< used to look for the local maximums in the energy vector */
	SignalPtr _filteredFirstDerivative; /**< first derivated filtered with a median filter */
	SignalPtr _secondDerivative; /**< used to look for the local maximums in the energy vector */
	std::vector<std::vector<unsigned int> > _microPairIdx; /**< contains the relationship between micro-pair idx and micros idx (k->{i.j}) */
	std::vector<std::shared_ptr<dsp::GeneralisedCrossCorrelation> > _gcc; /**< vector of GCC objects to compute the correlations (one object for each pair) */

	/**
	   * @brief allocate Allocates memory.
	   */
	void allocate();

	/**
	   * @brief generateLookupTable Generates a lookup table used to accelerate the real time processing.
	   */
	void generateLookupTable();

	/**
	   * @brief computeCorrelations  Computes the correlations for each micro pair using the GCC algorithm. Stores the
	   * result in _correlations.
	   * @param analysisFrames  Frames to be processed (FFT in CCS format).
	   * @param wienerCoefs  Wiener coeficients used in the GCC (not implemented yet).
	   */
	void computeCorrelations(const SignalVector &analysisFrames, SignalVector &wienerCoefs);

	/**
	   * @brief computeEnergyInDOA  Computes the energy in each DOA from the correlations computed previously.
	   * Stores the result in _energyInDOA.
	   */
	void computeEnergyInDOA();

	/**
	   * @brief selectDOA  Selects the numOfSources DOA's with maximum energy.
	   * @param DOA  Vector with the DOA estimations.
	   * @param prob  Probabilities associated with each DOA.
	   * @param numOfSources  Number of source to search.
	   */
	void selectDOA(SignalPtr DOA, SignalPtr prob, int numOfSources);
};


}

#endif
