/*
* Beamformer.h
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
#ifndef __MCA_BEAMFORMER_H_
#define __MCA_BEAMFORMER_H_

#include <mcarray/ArrayDescription.h>
#include <mcarray/mcadefs.h>

#include <memory>

namespace mca
{


class Beamformer
{

    public:

	Beamformer(int sampleRate, ArrayDescription microphonePositions, int fftCCSLength, unsigned int nchannels);
	virtual ~Beamformer(){}

	/**
	   * @brief processFrame  Processes one frame and generates the output frames combinating the input frames
	   * according to the specified DOA.
	   * @param inputAnalysisFrames  input frames to be processed (FFT in CCS format).
	   * @param outputFrame  output frame
	   * @param DOA  DOA to be pointed.
	   */
	void processFrame(SignalVector &inputAnalysisFrames, SignalPtr outputFrame, double DOA);

    private:

	int _sampleRate; /**< sample rate of the signals to be processed */
	int _fftCCSLength; /**< real length of the buffers containing the FFT in CCS format */
	int _nchannels; /**< number of channels in input signal (number of microphones) */
	ArrayDescription _microphonePositions; /**< description of the array (position of each microphone). */

	SignalPtr _phase; /**< phase ramp used to apply a delay to each channel */
	SignalCPtr _complexRamp; /**< complex "ramp" used to apply a delay to each channel */
	SignalPtr _ones; /**< used as magnitude to compute the complex ramp */
	SignalCPtr _channelSignal; /**< delayed signal of one channel, used to compute the output signal */

	/**
	   * @brief allocate Allocates memory.
	   */
	void allocate();
};


}

#endif
