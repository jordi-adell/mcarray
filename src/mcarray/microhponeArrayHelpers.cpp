/*
* microhponeArrayHelpers.cpp
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
#include <mcarray/microhponeArrayHelpers.h>
#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>

namespace mca{


// --- GCC class ------------------




// --- Helper functions ----------------------------


double getSpeedOfSound()
{
  // This is speed at 25ºC in standard atmosphere:
  // http://en.wikipedia.org/wiki/Speed_of_sound
  return 346.1;
}


float doaToDelayFarField(float doa, float microDist, bool useDegrees)
{
  //                                                * S
  //                                           .   ·
  //                            --   .            ·
  //                d       .                    ·
  //                .                           ·
  //180º  ------- (m)-----------D------------- (m) --------  0º
  //  t = d/c = (D/c)*cos(doa);
  //
  //  D = micro distance
  //  d = delay distance
  //  t = time distance (delay)
  //  doa = direction of arrival (-90 to 90) or (-pi/2 to pi/2)

  if (useDegrees)
  {
    doa = doa*M_PI/180;
  }
  float delay = (microDist * sin(doa)) / getSpeedOfSound();
  return delay;
}

float doaToDelayFarFieldSamples(float doa, float microDist, int sampleRate)
{
  return (doaToDelayFarField(doa, microDist)*sampleRate);
}

float delayToDOA(float delay, float microDist)
{
  return acos(delay*getSpeedOfSound()/microDist);
}

float delaySamplesToDOA(float delay, float microDist, float sampleRate)
{
  delay = delay/sampleRate;
  return delayToDOA(delay, microDist);
}

float maxFreqForSpatialAliasing(float microphoneDistance)
{
  //@TODO check the doubling of the distance
  return (getSpeedOfSound()/(2*microphoneDistance));
}

SignalPtr toDegrees(SignalPtr radians, int length)
{
  SignalPtr degrees;
  degrees.reset(new BaseType[length]);
  wipp::multC(180/M_PI, radians.get(), degrees.get(), length);

  return degrees;
}

SignalPtr toRadians(SignalPtr degrees, int length)
{
  SignalPtr radians;
  radians.reset(new BaseType[length]);
  wipp::multC(180/M_PI, radians.get(), degrees.get(), length);

  return radians;
}


float angle2DOAidx(float angle, float _doaStep)
{
  angle = std::max(static_cast<double>(angle), -M_PI_2);
  angle = std::min(static_cast<double>(angle),  M_PI_2);
  return (static_cast<int>((angle + M_PI_2)/_doaStep));
}

float doaIdx2angle(int idx, float _doaStep)
{
  return ((static_cast<float>(idx)*_doaStep)-M_PI_2);
}

double calculateBinauralPower(const BaseTypeC *left, const BaseTypeC *right, BaseTypeC *mixed, int length)
{
  BaseType magn[length];
  BaseType mean;

  const wipp::wipp_complex_t *w_left = reinterpret_cast<const wipp::wipp_complex_t*>(left);
  const wipp::wipp_complex_t *w_right = reinterpret_cast<const wipp::wipp_complex_t*>(right);
  wipp::wipp_complex_t *w_mixed = reinterpret_cast<wipp::wipp_complex_t*>(mixed);

  wipp::copyBuffer(w_left, w_mixed, length);
  wipp::add(w_right, w_mixed, length);
  wipp::magnitude(w_mixed, magn, length);
  wipp::sqr(magn, length);
  wipp::mean(magn, length, &mean);

  mean = 10*log10(mean/2);
  return mean;
}


}


