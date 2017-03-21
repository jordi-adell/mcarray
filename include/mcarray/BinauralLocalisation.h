/*
* BinauralLocalisation.h
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

#ifndef __BINAURALLOCALISATION_H
#define __BINAURALLOCALISATION_H

#include <mcarray/SoundLocalisationImpl.h>
#include <dspone/algorithm/gralCrossCorrelation.h>
#include <dspone/rt/ShortTimeFourierAnalysis.h>

namespace mca {


/** This class implements the Jefferss-Colburn model for binaural localisation
	 * as explained in
	 * Stern, Richard M, Evandro Gouvêa A, Chanwoo Kim, Kshitiz Kumar, and and Hyung-Min Park. 2008.
	 * BINAURAL AND MULTIPLE-MICROPHONE SIGNAL PROCESSING MOTIVATED BY AUDITORY PERCEPTION.
	 * In Proceedings of (HSCMA) Hands-Free Speech Communication and Microphone Arrays.
	 *
	 * Audio buffers are processed by calling the processBuffer function and the DOA result
	 * is set but calling the setDOA function of the callback supplied and derived from the
	 * abstract definition above.
	 **/
class TemporalGCCBinauralLocalisation : public SoundLocalisationImpl, public dsp::ShortTimeAnalysis
{

    public:
	TemporalGCCBinauralLocalisation(int sampleRate, ArrayDescription microphonePositions);

    private:

	// Interal parameter of the algorithm.
	static constexpr double _frameRate = 0.075; /**< the length of the frame to work with, in  seconds  */
	static constexpr double _doaMemoryFactor = 0.500; /**< the weight assign to prevous DOA in when updating _currentDOA */
	static constexpr double _doaMemoryFactorSilence = 0.800; /**< the weight assign to prevous DOA in when updating _currentDOA */

	const double _microphoneDistance; /**< distance between the two microphones used to record the audio signal, in meters */

	// Some values calculated in construction and that can not change over time.
	const int _sampleRate; /**< sample rate of the signals to be processed */
	const size_t _ndelays; /**< number of delays in the delay chain */
	const int _windowSize;
	const int _analysisLength;

	// Some internal variables
	double _signalPower;  /**< the power of the last frame used to calculated the DOA */

	// Some vectors used for internal operations.
	// Their memory is reserved on construction
	SignalPtr _leftDelayed; /**< left channel delayd by some samples */
	SignalPtr _rightDelayed; /**< right channel delayd by some samples */
	SignalPtr _crosCorr; /**< cross correlation of 2*_ndelay+1 samples between channels */
	SignalPtr _index; /**< vector of indicators, is an index of power in the DOA correspondin to each delay chain */
	SignalPtr _triangle; /**< a triangle of small magintude used to steer DOA to 0º in case of consecutive delays with same maximum index */
	SignalPtr _mixedChannel; /**< vector used to mix right and left channels */

	/**
	   * @brief frameAnalysis  implements the abstract class ShortTimeAnalysis unwindows each frame,
	   * and the analysis data is the frame itself.
	   * @param inFrame  input frame signal
	   * @param analysis  analysis of the frame. For this class the frame unwindowed version.
	   * @param frameLength  Length of the frame (in samples)
	   * @param analysisLength  Lenfth of the analysis (here the same as frameLength)
	   */
	virtual void frameAnalysis(BaseType *inFrame, BaseType *analysis, int frameLength, int analysisLength, int channel);

	/**
	   * @brief processParametrisation  Process the signal in _analysisFrame
	   * calculated by the frameAnalysis(...) function, caulcates the DOA and calls the callback
	   * with the obtained value
	   */
	virtual void processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
					    std::vector<double*> &dataChannels, int dataLength);

	/**
	   * @brief delay  Delays a frame by some samples
	   *     outFrame[m] = inFrame[m-delay]
	   * @param inFrame  frame to be delayed
	   * @param inLength  length of the input frame
	   * @param outFrame  vector to store the delayed frame in
	   * @param outLength  length of the delays frame.
	   * @param delay  number of samples to delay the frame.
	   */
	void delay(BaseType *inFrame, int inLength, BaseType *outFrame, int outLength, int delay);

	/**
	   * @brief crossCorrelation  calculates the cross-correlation between two signal frames.
	   * The correlation is calculate for lags in the range [-outLength/ outLength/2]
	   * @param leftFrame  pointer to the left channel frame.
	   * @param leftLength  lengthe of the left channel frame.
	   * @param rightFrame  pointer to the right channel frame
	   * @param rightLength  length of the right channel frame.
	   * @param outFrame  pointer to the output frame where corss-correlation is stored.
	   * @param outLength  length of the output frame
	   */
	void crossCorrelation(BaseType *leftFrame, int leftLength, BaseType *rightFrame, int rightLength, BaseType *outFrame, int outLength);

	/**
	   * @brief correlation  Calculates correlation between to frames by simply multiplying and
	   * add the input frames:
	   *   c = (1/N)*sum_i (l[i]*r[i])
	   * @param leftFrame  pointer to the left channel frame. (l)
	   * @param rightFrame pointer to the right channel frame. (r)
	   * @param length  length of both frames  (N)
	   * @return  returns the correlation between input frames (c)
	   */
	double correlation(BaseType *leftFrame, BaseType *rightFrame, int length);
	/**
	   * @brief integration  integrates a vector (sum all its componenets).
	   * This function is used to integrate the crossCorrelation of both channels.
	   * @param inFrame  frame to integrate
	   * @param length  length of the input frame.
	   * @return   the result of the integration.
	   */
	double integration(BaseType *inFrame, int length);
	/**
	   * @brief maxminNormalisation  nomalised a frame by add its minimum to it
	   * and dividing it by its maximum after the previous addition. The result
	   * is stored in the input vector itself.
	   * @param vector  pointer to the vector to be normalised.
	   * @param length  length of the input vector.
	   */
	void maxminNormalisation(BaseType *vector, int length);

	/**
	   * @brief estimateDOALogLikelihood  this function estimates the log-likelihood of the
	   * the DOA stored in _index[maxindex] and stores the value in prob:
	   *     log(_index[maxindex] / sum_i abs(_index[i]))
	   * This is not a good estimatino, some research need to be to improve this
	   * estimation.
	   * @param maxindex  index of the delay in _index corresponding to the DOA
	   * we want to calculate the log-likelyhood.
	   * @param prob  The log-likelihood is stored in this variable.
	   */
	void estimateDOALogLikelihood(int maxindex, BaseType &prob);

	/**
	   * @brief degrees2Samples  converts the DOA in degrees to the corresponding
	   * delay in samples.
	   * @param degrees  degres of the DOA
	   * @return returns the corresponding delay in samples.
	   */
	int degrees2Samples(double degrees, int sampleRate);
	/**
	   * @brief samples2Degrees  converts the observed delay in samples between microphones
	   * to degrees of the DOA
	   * @param delay   number of samples indicating the delay between both microphones.
	   * @return  returns de DOA in degrees.
	   */
	double samples2Degrees(int delay, int ndelays);

	/**
	   * @brief setPowerFloor  Sets the power floor and when the necessary amount of samples have been
	   * consumed set the variable __noiseEstimated to true and the variable _powerFloor.
	   * Used calculatePower to obtain the power value.
	   * to the estimated value.
	   * @param analysisFrames  input frames
	   * @param length  length of the input vectors
	   * @param nchannels  number of channels in analysisFrames
	   * @param sampleRate  sampleRate
	   * @return   returns the power calculated.
	   */
	virtual BaseType setPowerFloor(SignalVector &analysisFrames, int analysisLength, int nchannels, int sampleRate);
	BaseType setPowerFloor(std::vector<BaseType*> &analysisFrames, int analysisLength, int nchannels, int sampleRate);

};


class FreqGCCBinauralLocalisation : public SoundLocalisationImpl, public dsp::STFTAnalysis
{
    public:
	FreqGCCBinauralLocalisation(int sampleRate, ArrayDescription microphonePositions, bool usePowerFloor=true);
	virtual void setProbability(const double *doas, double *probs, int size);

    private:

	static constexpr float _frameRate = 0.075; /**< related with the length of the frame to work with, in  seconds  */
	static constexpr float _noiseMarginDB = 6.0; /**< def: 3 margin over the noise level to decide that signal is present  */
	static constexpr float _maxCorrMemoryFactor = 0.8; /**< max memory factor to smooth the correlation (weigth assigned to the previous correlation).  */
	static constexpr float _maxDoaMemoryFactor = 0.6; /**< max memory factor used to smooth the DOA values (weigth assigned to the previous DOA).   */
	float _corrMemoryFactor; /**< Current memory factor used to smooth the correlation. It goes to zero in silence periods. */
	float _doaMemoryFactor; /**< Current memory factor used to smooth the DOA values. It goes to zero in silence periods.  */
	const double _microphoneDistance; /**< distance between the two microphones used to record the audio signal, in meters */
	int _silenceFramesCounter; /**< Counter of the consecutive frames in silence (without signal). Used to control memory factors.  */
	int _sampleRate; /**< sample rate of the signals to be processed */
	const float _doaStep; /**< Step in degrees between DOA values used in the DOA grid. It determines the resolution.  */
	const int _numSteps; /**< Number of steps in DOA grid (calculated from doaStep).  */
	bool _usePowerFloor; /**< Flag that indicates if the power floor is used to detect the presence of signal.  */
	std::unique_ptr<uint16_t> _sampledDOAs;
	SignalPtr _magnitude; /**< vector used to compute the magnitude of the signal. */
	SignalPtr _power; /**< vector used to compute the power of the signal. */
	SignalCPtr _correlations; /**< vector used to store the correlation. */
	SignalPtr _correlationsReal; /**< vector used to store the real part of the correlation. */
	SignalPtr _prevCorrelationsReal; /**< vector used to store the previous correlation. */
	SignalPtr _triangle; /**< triangle used to give more weigth to the center DOAs in the correlation. */
	SignalCPtr _mixedChannel; /**< vector used to mix right and left channels */
	SignalPtr _samplesDelay; /**< Delay (in samples) used to compute the correlation, according to doaStep and numSteps.  */

	dsp::GeneralisedCrossCorrelation _gcc; /**< pointer to the object which implements the GCC algorithm. */

	/**
	   * @brief processParametrisation  Process the signal in _analysisFrames
	   * calculated by the frameAnalysis(...) function, caulcates the DOA and calls
	   * the callback with the obtained value.
	   */
	virtual void processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
					    std::vector<double*> &dataChannels, int dataLength);
	/**
	   * @brief setPowerFloor  Sets the power floor and when the necessary amount of samples have been
	   * consumed set the variable __noiseEstimated to true and the variable _powerFloor.
	   * Used calculatePower to obtain the power value.
	   * to the estimated value.
	   * @param analysisFrames  input frames
	   * @param length  length of the input vectors
	   * @param nchannels  number of channels in analysisFrames
	   * @param sampleRate  sampleRate
	   * @return   returns the power calculated.
	   */
	virtual BaseType setPowerFloor(SignalVector &analysisFrames, int analysisLength, int nchannels, int sampleRate);
	BaseType setPowerFloor(std::vector<BaseType*> &analysisFrames, int analysisLength, int nchannels, int sampleRate);


	/**
	   * @brief For debuggin purposes
	   */
	void debugPlot(int idx);

};






}


#endif // __BINAURALLOCALISATION_H
