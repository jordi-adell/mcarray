/*
* mcalogger.h
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

#ifndef __DSPONE_LOGGER_H_
#define __DSPONE_LOGGER_H_

#include <iostream>
#include <vector>


#ifdef MCA_LOGGER

#include <mcarray/mcaloggerImpl.h>

#else

#ifndef INFO_STREAM
#define INFO_STREAM(x)
#endif

#ifndef WARN_STREAM
#define WARN_STREAM(x)
#endif

#ifndef ERROR_STREAM
#define ERROR_STREAM(x)
#endif

#ifndef DEBUG_STREAM
#define DEBUG_STREAM(x)
#endif

#ifndef TRACE_STREAM
#define TRACE_STREAM(x)
#endif

#ifndef ERROR_STREAM_ONCE
#define ERROR_STREAM_ONCE(x)
#endif

#ifndef WARN_STREAM_ONCE
#define WARN_STREAM_ONCE(x)
#endif

#ifndef INFO_STREAM_ONCE
#define INFO_STREAM_ONCE(x)
#endif

#ifndef DEBUG_STREAM_ONCE
#define DEBUG_STREAM_ONCE(x)
#endif

#ifndef SET_LOG_LEVEL
#define SET_LOG_LEVEL(x)
#endif

#endif


#endif //__DSP_LOGGER_H_
