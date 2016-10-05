/*
* mcalogger.cpp
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
#define _LOGGER

#include <mcarray/mcalogger.h>
#include <boost/filesystem.hpp>

namespace  mca {

	Logger::Logger(LogLevel level, bool enable, std::ostream &os) : _level(level), _os(os), _enabled(enable)
	{
	  //	  \x1b[32m
	    _levelNames = {"FATAL", "ERROR", "WARN ", "INFO ", "DEBUG", "TRACE"};
	}

	std::ostream &Logger::log(LogLevel level, const char *file, int line)
	{
	    std::string fname = boost::filesystem::basename(file);
	    if (_enabled && level <= _level)
	    {
		_os << "[" << _levelNames.at(level) << "] " <<  "[" << fname << ":" << line <<  "] ";
		return _os;
	    }
	}

	bool Logger::isEnabled(LogLevel level) { return (_enabled && level <= _level); }
	void Logger::enable() { _enabled = true;}
	void Logger::disable() { _enabled = false; }

}

