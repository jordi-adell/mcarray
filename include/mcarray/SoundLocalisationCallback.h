/*
* SoundLocalisationCallback.h
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
#ifndef SOUNDLOCALISATIONCALLBACK_H
#define SOUNDLOCALISATIONCALLBACK_H

#include <mcarray/mcadefs.h>
#include <dspone/dsp.h>

namespace mca {


/** This abstract class is used as a callback to recieve the DOA values
	 * at the time they are calculated along with some suplementary
	 * information such as the probability and the singal power.
	 * If you want to collect such data  you need to implement
	 * a derived class from this one.
	 **/
class LocalisationCallback
{
    public:
	LocalisationCallback(){}
	virtual ~LocalisationCallback(){}
	/**
	   * @brief setDOA This function is called bt the BinauralLocalisation
	   * class every time a new doa is calculated.
	   * @param doa  Direction Of Arrival (DOA) of the detected sources in degrees.
	   * The range of this value is [-90º 90º]
	   * 0º indicate front direction, 90º indicates left and -90º indicates right.
	   * @param prob  The the probability of the observated DOAs.
	   * @param power  The power of the frame used to calculate the DOAs.
	   * If this power is lower that a floor power calculated at start-up, this funciton is
	   * not called.
	  */
	virtual void setDOA(SignalPtr doa, SignalPtr prob, double power, int numOfSources) = 0;
};


/**
	 * @brief The DummyLocalisationCallback class  a Dummy implementation of a
	 * Localisation Callback. It just prints out the set DOA.
	 */
class DummyLocalisationCallback : public LocalisationCallback
{
    public:
	DummyLocalisationCallback(){}
	~DummyLocalisationCallback(){}
	/**
	   * @brief setDOA  fucntion to called by a BinauralLocalisation object.
	   * @param doa  Direction Of Arrival.
	   * @param prob  probavility assigned to the obtained DOA.
	   * @param power  power of the frame used to calculated the DOA.
	   */
	using LocalisationCallback::setDOA;
	virtual void setDOA(double doa, double prob, double power);
};

}


#endif // SOUNDLOCALISATIONCALLBACK_H
