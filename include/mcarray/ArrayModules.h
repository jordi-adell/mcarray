/*
* ArrayModules.h
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
#ifndef ARRAYMODULES_H
#define ARRAYMODULES_H

#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/ArrayDescription.h>
#include <dspone/dsp.h>
#include <dspone/rt/ShortTimeAnalysis.h>

#include <memory>

namespace mca {

class LocalisationCallback;

/**
       * @brief The SoundLocalisation class
       * This class computes the Direction of Arrival (DOA) of a sound source.
       * And Publishes it. It has not output ports.
       */
class SoundLocalisation : public dsp::SignalAnalyser
{
    public:
	/**
	 * @brief SoundLocalisation
	 * @param sampleRate   sample rate of the input signal.
	 * @param microphonePositions   Description of the microphone array.
	 * @param callback    Callback to call, whenever a new DOA is computed.
	 */
	SoundLocalisation(int sampleRate,
			  ArrayDescription microphonePositions,
			  LocalisationCallback *callback=NULL);
	virtual ~SoundLocalisation();
    private:
	std::unique_ptr<dsp::ShortTimeAnalysis> _impl;
};



/**
       * @brief The NoiseReductionEM class applies an algorithm for noise reduction
       * bases on PSD noise estimation and a Wiener filter based on SNR.
       * Noise is estimated and reduced for each channel independently.
       */
class BinauralMasking: public dsp::SignalProcessor
{
    public:
	/**
	   * This type is for identifying masking methods
	   * FACTOR = divides sptatial masked time-frequency bins by a factor,
	   *          and temporal masked by another factor
	   * RELATIVE = Uses the low pass filtered power of the time-frequency bin
	   *            with respect to the current power to calculate the factor
	   *            to multiply the rejected signal.
	   * FULL  =  Sets the rejected time-frequency bin to zero (removes it completlly)
	   * NOISY = Calculate the factor to multiply the rejected signal to obtain
	   *         a flat spectrum in order to eliminated it later with the noise
	   *         reduction system
	   **/

	typedef enum {FACTOR =0, RELATIVE=1, FULL=3, NOISY=4, NOTHING=5} MaskingMethod;

	/**
	    * This type identifies which part of the algorithm is actually exectued.
	    * SPATIAL = frames will only be masked following spatial criteria.
	    * TEMPORAL = frmaes will only be masked following temportal criteria (to reduce reverberation)
	    * BOTH = both previous criteria are followed.
	    **/
	typedef enum {BOTH=0, SPATIAL=1, TEMPORAL=2} MaskingAlg;

	/**
	 * @brief BinauralMasking
	 * @param sampleRate  sample rate of the signal to be processed.
	 * @param nchannels  number of channels to be process.
	 */
	BinauralMasking(int samplerate,
			ArrayDescription microphones,
			float lowFreq = 400,
			float highFreq = 4000,
			MaskingMethod mmethod = RELATIVE,
			MaskingAlg algorithm = BOTH);

	virtual ~BinauralMasking();
    private:
	std::unique_ptr<dsp::ShortTimeProcess> _impl;
};



}



#endif // ARRAYMODULES_H
