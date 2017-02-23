/*
* mcaloggerImpl.h
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
#ifndef __MCA_LOGGER_IMPL_H_
#define __MCA_LOGGER_IMPL_H_

#include <iostream>
#include <vector>
#include <memory>

#ifdef MCA_LOGGER
#if MCA_LOGGER
#undef MCA_LOGGER
#define MCA_LOGGER FATAL
#endif
#endif

namespace mca
{

class Logger
{
    public:
	typedef enum{FATAL  = 0, ERROR = 1, WARNING = 2, INFO = 3, DEBUG = 4, TRACE = 5} LogLevel;

	Logger(LogLevel level = FATAL, bool enable = true, std::ostream &os = std::cout);
	static Logger& logger(LogLevel level = FATAL, bool enable = true, std::ostream &os = std::cout);

	std::ostream &log(LogLevel level, const char *file, int line);

	bool isEnabled(LogLevel level);
	void enable();
	void disable();

	void setLevel(std::string level);
	void setLevel(LogLevel level);

    private:
	LogLevel _level;
	std::ostream &_os;
	bool _enabled;
	std::vector<std::string> _levelNames;

	static std::unique_ptr<Logger> _logger;
};


}

#define _logger_ mca::Logger::logger(mca::Logger::MCA_LOGGER)

#define LOG_STREAM(level, what) {if (_logger_.isEnabled(level)) _logger_.log(level, __FILE__, __LINE__) << what << std::endl;}
#define DEBUG_STREAM(what) LOG_STREAM(mca::Logger::DEBUG,   what)
#define ERROR_STREAM(what) LOG_STREAM(mca::Logger::ERROR,   what)
#define WARN_STREAM(what)  LOG_STREAM(mca::Logger::WARNING, what)
#define INFO_STREAM(what)  LOG_STREAM(mca::Logger::INFO,    what)
#define TRACE_STREAM(what) LOG_STREAM(mca::Logger::TRACE,   what)

#define WARN_STREAM_ONCE(what)  {static bool printed=false; if (!printed) LOG_STREAM(mca::Logger::WARNING, what); printed=true;}
#define ERROR_STREAM_ONCE(what) {static bool printed=false; if (!printed) LOG_STREAM(mca::Logger::ERROR, what);   printed=true;}
#define INFO_STREAM_ONCE(what)  {static bool printed=false; if (!printed) LOG_STREAM(mca::Logger::INFO, what);    printed=true;}
#define DEBUG_STREAM_ONCE(what) {static bool printed=false; if (!printed) LOG_STREAM(mca::Logger::DEBUG, what);   printed=true;}
#define TRACE_STREAM_ONCE(what) {static bool printed=false; if (!printed) LOG_STREAM(mca::Logger::TRACE, what);   printed=true;}

#define SET_LOG_LEVEL(level) _logger_.setLevel(level)

#endif //__MCA_LOGGER_H_
