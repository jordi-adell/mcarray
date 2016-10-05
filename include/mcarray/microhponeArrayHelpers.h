/*
* microhponeArrayHelpers.h
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
#ifndef __MICROHPONEARRAYHELPERS_H
#define __MICROHPONEARRAYHELPERS_H

#include <mcarray/mcadefs.h>
#include <dspone/complex.h>

#include <boost/shared_array.hpp>

#include <vector>
#include <math.h>


namespace mca
{


double getSpeedOfSound();
float doaToDelayFarField(float doa, float microDist, bool useDegrees=false);
float doaToDelayFarFieldSamples(float doa, float microDist, int sampleRate);
float delayToDOA(float delay, float microDist);
float delaySamplesToDOA(float delay, float microDist, float sampleRate);
float maxFreqForSpatialAliasing(float microphoneDistance);
inline float toDegrees(float radians){return radians*M_1_PI*180;}
SignalPtr toDegrees(SignalPtr radians, int length);
inline float toRadians(float degrees){return degrees*M_PI/180;}
SignalPtr toRadiasn(SignalPtr degrees, int length);
float angle2DOAidx(float angle, float _doaStep);
float doaIdx2angle(int idx, float _doaStep);

double calculateBinauralPower(const Complex *left, const Complex *right, Complex *mixed, int length);

}

#endif // __MICROHPONEARRAYHELPERS_H
