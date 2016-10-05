/*
* SoundLocalisationCallback.cpp
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
#include <mcarray/SoundLocalisationCallback.h>
#include <mcarray/mcalogger.h>

#include <sstream>


namespace mca {

// -------- DummyCallback -----------------------------------------------

void DummyLocalisationCallback::setDOA(double doa, double prob, double power)
{
    // The results are printed for debugging purposes.
    std::ostringstream oss;
    double degreeStep = 10;
    for (double angle = -90; angle < 90; angle += degreeStep)
    {
	if (angle - degreeStep/2 < doa && doa < angle + degreeStep/2)
	    oss << "*";
	else if (-degreeStep/2 < angle && angle < degreeStep/2)
	    oss << "|";
	else
	    oss << "-";
    }

    oss << " (" << doa << "º)";

    if (doa < -20)
	oss << "\tRIGHT";
    else if (doa > 20)
	oss << "\tLEFT";
    else
	oss << "\tFRONT";

    oss << " (p: "<< prob << ", P: " << power << ")";
    DEBUG_STREAM(oss.str());
    oss.str("");
}

}

