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

#include <mcarray/mcaloggerImpl.h>
#include <boost/filesystem.hpp>


namespace  mca {

std::unique_ptr<Logger> Logger::_logger;

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


void Logger::setLevel(LogLevel level)
{
  _level = level;
}

void Logger::setLevel(std::string level)
{
  auto it = std::find(_levelNames.begin(), _levelNames.end(), level);
  if (it != _levelNames.end())
  {
    setLevel(static_cast<LogLevel>(it - _levelNames.begin()));
  }
}

bool Logger::isEnabled(LogLevel level) { return (_enabled && level <= _level); }
void Logger::enable() { _enabled = true;}
void Logger::disable() { _enabled = false; }

Logger &Logger::logger(LogLevel level, bool enable, std::ostream &os)
{
  if (!_logger)
  {
    _logger.reset(new Logger(level, enable, os));
  }
  return *(_logger.get());
}

}

