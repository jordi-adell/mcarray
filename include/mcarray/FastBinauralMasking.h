/*
* FastBinauralMasking.h
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
#ifndef __FAST_BINAURALMASKING_H
#define __FAST_BINAURALMASKING_H

#include <mcarray/ArrayModules.h>
#include <dspone/algorithm/gralCrossCorrelation.h>
#include <dspone/rt/ShortTimeFourierTransform.h>
#include <dspone/filter/FilterBank.h>

#include <boost/scoped_array.hpp>


namespace mca {


/** Binaural spatial-temporal masking as can be found in:
	*
	* Kim, C., K. Kumar, and R. M. Stern. 2011.
	* "Binaural sound source separation motivated by auditory processing."
	* In IEEE International Conference on Acoustics, Speech, and Signal Processing, 4.
	*
	* There are some differences in the algorithm is that we enhance accepted time-frequency bins,
	* and not only degrade the regected ones. This is because we are working with low
	* power signals.
	*
	* Also we implemente thre masking functions a part from the one proposed in the paper,
	* here we can choose to completly remove the rejected bin, or to mutiply the rejected
	* signal by factor for spatial masking and by a different one for temporal masking.
	* This can be useful in cases were we rely more on one masking type.
	*
	* Use use a mel-scaled filter bank instead of a gammtone filterbank.
	*
	**/

class FastBinauralMasking : public dsp::STFT
{
    public:

	typedef BinauralMasking::MaskingMethod MaskingMethod;
	typedef BinauralMasking::MaskingAlg MaskingAlg;

	/**
	   * @brief FastBinauralMasking constructor
	   * @param nchannels  number of channels to be process
	   * (currently Binaural Masking only accepts two input channels, otherwise
	   * it would not be binaural! )
	   * @param samplerate
	   * @param microDistance  distance between the microphones in metres
	   * @param mmethod  set the masking method used
	   */
	FastBinauralMasking(int samplerate,
			    double microDistance,
			    float lowFreq,
			    float highFreq,
			    MaskingMethod mmethod = BinauralMasking::RELATIVE,
			    MaskingAlg algorithm = BinauralMasking::BOTH);

	/**
	   * @brief Processes the analysis buffers and changes the time-frequency bins stored in the analysis buffer
	   * accoording to the implemented algorithm. Performes the masking itself.
	   */
	virtual void processParametrisation(std::vector<double*> &analysisFrames, int analysisLength,
					    std::vector<double*> &dataChannels, int dataLength);

	/**
	   * @brief getNonMaskingAngle
	   * @return the angle in degrees for the region where signal is accepted.
	   */
	inline int getNonMaskingAngle(){return _phi*180*M_1_PI;}
	/**
	   * @brief getMicroPhoneDistance
	   * @return returns the distance between microphones in metres.
	   */
	inline float getMicroPhoneDistance(){return _microDistance;}
	/**
	   * @brief getSpatialMaskingFactor
	   * @return  the factor used for spatial masking in the FACTOR masking method
	   */
	inline float getSpatialMaskingFactor(){return 1/_spatialMaskingFactor;}
	/**
	   * @brief getTemporalMaskingFactor
	   * @return the factor used for temporal masking in the FACTOR masking method
	   */
	inline float getTemporalMaskingFactor(){return 1/_temporalMaskingFactor;}

    private:

	/**
	   * @brief Parameters of the algorithm
	   */
	static constexpr  int _nBins = 45;
	static constexpr float _frameRate = 0.050; // in secods, will be windowShift and half of windowSize
	static constexpr double _phi = 10*M_PI/180; // Degrees to radians
	static constexpr float _forgetingFactor = 0.04; // necessary for memory in temporal masking

	// Parameters for FACTOR method
	static constexpr float _temporalMaskingFactor=3;  // Masked signal is divided by this factor (3 ~ -10dB)
	static constexpr float _spatialMaskingFactor=10;  // Masked signal is divided by this factor (10 ~ -20dB)

	// Parameter for enhanced frames
	static constexpr float _enhanceFactor = 1; // enhanced time-frequency bins are multiplied by this factor (2 ~ 6dB)

	// Parameters for RELATIVE method
	static constexpr float _scalingFactor = 0.01; // ro in Kim et al.'s paper (they say 0.01 is a good value) (0.01 ~ -40dB)

	// Parameter to scale the low-passed filtered input power to decide if it is masked by temporal masking. With a factor of 1 the signal
	// is attenuated before the arribal of the reverberation. With 0.999 we observed a large improvement.
	static constexpr float _rejectTemporalFactor = 0.999;


	/**
	   * Configuration parameters
	   **/
	const int _sampleRate;
	const double _microDistance;
	const MaskingMethod _mmethod;
	const MaskingAlg _algorithm;
	// Pre filter parameters
	const float _minFreq;
	const float _maxFreq;
	const int _nchannels;
	const int _fftOrder;
	const int _oneSidedFFTLength;
	const int _windowSize;

	int _firstCall;

	//          boost::scoped_ptr<GeneralisedCrossCorrelation> _gcc;


	/**
	   * @brief Mel-scaled filter bank.
	   */
	std::unique_ptr<dsp::FilterBank> _filterBank;
	int _coeficientsLength;
	boost::scoped_array<BaseTypeC> _filterCoeficients;
	boost::scoped_array<BaseType> _fftLeftFrame;
	boost::scoped_array<BaseType> _fftRightFrame;
	boost::scoped_array<BaseType> _outLeftFrame;
	boost::scoped_array<BaseType> _outRightFrame;



	/**
	   * @brief Low pass-filtered power.
	   * Memory of the temporal masking.
	   */
	boost::scoped_array<BaseType> _shortTimePower; // A value for each time-frequency bin

	/**
	   * @brief Low pass-filtered power.
	   * Memory of the temporal masking.
	   */
	boost::scoped_array<BaseType> _noiseEstimatePower; // A value for each time-frequency bin

	/**
	   * @brief Normalised correlation thresholds for accptance/rejection
	   * in spatial masking.
	   */
	std::vector<double> _thresholds;

	/**
	   * Deprecated
	   ***/
	double localise(BaseType* left, BaseType *right, int length);

	/**
	   * @brief Calculates the normalised correlation used in spatial masking.
	   * @param left   left channel buffer
	   * @param right  right channel buffer
	   * @param length  number of smaples in each channel
	   * @return   the nromalised correlation.
	   */
	double normaliseCorrelation(BaseType* left, BaseType *right, int length);
	double normaliseFFTCorrelation(BaseType* left, BaseType *right, int length);
	double generalisedCrossCorrelation(BaseType* left, BaseType *right, int length);

	/**
	   * @brief Decides whethe the time-frequency bin has to be masked or not
	   * acording to spatial information.
	   * @param left  left channel signal corresponding to the filterd singal in
	   * the specified bin.
	   * @param right  right channel signal corresponding to the filterd singal in
	   * the specified bin.
	   * @param length  length of left and right buffers
	   * @param bin  bin to which the given signal corresponds.
	   * @return  true if the bin has to be masked.
	   */
	inline bool spatialMasking(BaseType *left, BaseType *right, int length, int bin);

	/**
	   * @brief Decides whethe the time-frequency bin has to be masked or not
	   * acording to temporal information.
	   * @param left  left channel signal corresponding to the filterd singal in
	   * the specified bin.
	   * @param right  right channel signal corresponding to the filterd singal in
	   * the specified bin.
	   * @param length  length of left and right buffers
	   * @param bin  bin to which the given signal corresponds.
	   * @return  true if the bin has to be masked.
	   */
	inline bool temportalMasking(BaseType *left, BaseType *right, int length, int bin);

	/**
	   * @brief Calculate the power (linear scale) for the specified frame
	   * @param left  left channel buffer
	   * @param right  right channel buffer
	   * @param length  length of the buffers
	   * @return the power of the averaged channels
	   */
	inline BaseType getFramePower(BaseType *left, BaseType *right, int length);

	/**
	   * @brief getFRamePower  same as previous function but for a single channel
	   * @param frame  signal to get power
	   * @param length  length of the buffer
	   * @return the power of the channel
	   */
	inline BaseType getPower(BaseType *frame, int length);

	/**
	   * @brief Executed a masking method over the given frame
	   * @param frame  signal buffer
	   * @param length  length of frame
	   * @param factor  factor to be applied to the frame. This paramters is only used
	   * if the specified method used a given factor.
	   * @param bin  bin at which the given signal belongs to.
	   */
	inline void maskFrame(BaseType *frame, int length, float factor, int bin);


	/**
	   * @brief Enhance a frame in the case it is not rejected by
	   * multiplying it by a factor (has to be greater than one)
	   * @param frame  signal buffer
	   * @param length   length of the buffer
	   */
	inline void enhanceFrame(BaseType *frame, int length);

	/**
	   * @brief Sets the frame buffer to zero.
	   * @param frame frame buffer
	   * @param length   length of the buffer
	   */
	void zeroFrame(BaseType *frame, int length);

	/**
	   * @brief Sets the frame buffer to ******************
	   * @param frame frame buffer
	   * @param length   length of the buffer
	   * @param bin  bin at which the given signal belongs to
	   */
	void noisyFrame(BaseType *frame, int length, int bin);

	/**
	   * @brief Multiplies the frame by a factor that is calculated
	   * as reported by Kim et al. (2011).
	   * @param frame buffer frame
	   * @param length  legth of the buffer
	   * @param bin  time-frequency bin the signal belongs to.
	   */
	void maskFrameByScaling(BaseType *frame, int length, int bin);

	/**
	   * @brief Divides the signal frame by the supplied factor
	   * @param frame  frame signal buffer
	   * @param length  length of the buffer
	   * @param factor  factor by whish the signal is divided
	   */
	void maskFrameByFactor(BaseType *frame, int length, float factor);

	/**
	   * @brief Calculates the normalised correlation thresholds for spatial
	   * masking and stores them in _thresholds
	   */
	void calculateThresholds();

	/**
	   * @brief Initialises variables and allocates memory.
	   */
	void init();


};

}


#endif // __FAST_BINAURALMASKING_H
