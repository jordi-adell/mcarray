/*
* ArrayDescription.cpp
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

#include <mcarray/ArrayDescription.h>
#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/mcarray_exception.h>

#include <mcarray/mcalogger.h>

#include <math.h>

namespace mca {

ArrayDescription::ArrayDescription()
{
}

ArrayDescription::~ArrayDescription()
{
}

ArrayDescription ArrayDescription::make_linear_array_description(const std::vector<double> &x)
{
  ArrayDescription desc;
  for (size_t i = 0; i < x.size(); ++i)
  {
    desc.pushPosition(x.at(i), 0, 0);
  }
  return desc;
}


double ArrayDescription::distance(const std::string &i, const std::string &j) const
{
  return distance(getId(i), getId(j));
}

double ArrayDescription::distance(ElementId i, ElementId j) const
{
  return sqrt(
	pow(getX(j) - getX(i),2) +
	pow(getY(j) - getY(i),2) +
	pow(getZ(j) - getZ(i),2)
	);
}

double ArrayDescription::distance(ElementList::const_iterator &i,
				  ElementList::const_iterator &j) const
{
  return sqrt(
	pow(j->second.get<X>() - i->second.get<X>(),2) +
	pow(j->second.get<Y>() - i->second.get<Y>(),2) +
	pow(j->second.get<Z>() - i->second.get<Z>(),2)
	);
}


double ArrayDescription::maxDistance() const
{
  ElementList::const_iterator itI, itJ;
  double maxDistance = 0;
  for (itI = _list.begin(); itI != _list.end(); ++itI)
  {
    for (itJ = _list.begin(); itJ != _list.end(); ++itJ)
    {
      if (itI == itJ) continue;
      double d = distance(itI, itJ);
      maxDistance = std::max(d, maxDistance);
    }
  }
  return maxDistance;
}

double ArrayDescription::minDistance() const
{
  ElementList::const_iterator itI, itJ;
  double minDistance = 0;
  for (itI = _list.begin(); itI != _list.end(); ++itI)
  {
    for (itJ = _list.begin(); itJ != _list.end(); ++itJ)
    {
      if (itI == itJ) continue;
      double d = distance(itI, itJ);
      minDistance = std::min(d, minDistance);
    }
  }
  return minDistance;
}


ArrayDescription::ElementId ArrayDescription::pushPosition(const ArrayPosition &position,
							   const std::string &name)
{
  ElementId id = getNewId();
  std::string name_ = name;
  if (name_.empty())
    name_ = std::to_string(id);

  if (_names.find(name_) != _names.end())
    throw(MCArrayException("Array element name aready used"));
  else
  {
    _names[name_] = id;
    _list[id] = position;
  }
  return id;
}

ArrayDescription::ElementId ArrayDescription::pushPosition(double x, double y, double z,
							   const std::string &name)
{
  ArrayPosition position(x, y, z);
  return pushPosition(position, name);
}


size_t ArrayDescription::size() const
{
  return _list.size();
}

void ArrayDescription::getPosition(const std::string &name, ArrayPosition &position) const
{
  getPosition(getId(name), position);
}

void ArrayDescription::getPosition(ElementId id, ArrayPosition &position) const
{
  if (_list.find(id) != _list.end())
    position = _list.at(id);
  else
    ERROR_STREAM("No id " << id << " found in the array list");
}

double ArrayDescription::getX(const std::string &name) const
{
  return getX(getId(name));
}

double ArrayDescription::getY(const std::string &name) const
{
  return getY(getId(name));
}

double ArrayDescription::getZ(const std::string &name) const
{
  return getZ(getId(name));
}

double ArrayDescription::getX(const ElementId &id) const
{
  return get(id, X);
}

double ArrayDescription::getY(const ElementId &id) const
{
  return get(id, Y);
}

double ArrayDescription::getZ(const ElementId &id) const
{
  return get(id, Z);
}

double ArrayDescription::get(const ElementId &id, AxisName axis) const
{
  if (_list.find(id) == _list.end())
    ERROR_STREAM("No element with id " << id << " present in the array description");
  const ArrayPosition &position = _list.at(id);

  double value = 0;
  switch(axis)
  {
    case X:
      value = position.get<X>();
    break;
    case Y:
      value = position.get<Y>();
    break;
    case Z:
      value = position.get<Z>();
  }
  return value;
}

void ArrayDescription::getX(std::vector<double> &x) const
{
  x.clear();
  for (ElementId i = 0; i < size(); ++i)
  {
    x.push_back(getX(i));
  }
}

void ArrayDescription::getY(std::vector<double> &y) const
{
  y.clear();
  for (ElementId i = 0; i < size(); ++i)
  {
    y.push_back(getY(i));
  }
}

void ArrayDescription::getZ(std::vector<double> &z) const
{
  z.clear();
  for (ElementId i = 0; i < size(); ++i)
  {
    z.push_back(getZ(i));
  }
}

std::string ArrayDescription::getName(ElementId id) const
{
  ElementNames::const_iterator it;
  for (it = _names.begin(); it != _names.end(); ++it)
  {
    if (it->second == id)
    {
      return it->first;
    }
  }
  ERROR_STREAM("Id not found in the array description names list");
  return "";
}

ArrayDescription::ElementId ArrayDescription::getId(const std::string &name) const
{
  ElementId id = -1;
  if (_names.find(name) != _names.end())
    id = _names.at(name);
  else
    ERROR_STREAM("No position with name " << name << " found in the array");
  return id;
}

ArrayDescription::ElementId ArrayDescription::getId(const ElementId &id) const
{
  ElementId id_ = -1;
  if (_list.find(id) == _list.end())
    ERROR_STREAM("No element with id " << id << " present in the array description");
  else
    id_ = id;
  return id_;
}

ArrayDescription::ElementId ArrayDescription::getNewId()
{
  ElementId id = -1;
  for (id = size(); id < size() + 1e10 && id > 0;++id)
  {
    if (_list.find(id) == _list.end())
    {
      break;
    }
  }
  if (id == -1)
  {
    throw(MCArrayException("Unable to find a free array id"));
  }
  return id;
}


bool ArrayDescription::empty() const
{
  return _list.empty();
}


double ArrayDescription::getBandwidth() const
{
  double maxD = maxDistance();
  if (maxD <= 0)
    return 0;
  else
    return (getSpeedOfSound()/(2*maxDistance()));
}


std::ostream& operator<<(std::ostream& os, const ArrayDescription& desc)
{
  for (int i = 0; i < desc.size(); ++i)
  {
    os << "{"<<  desc.getName(i) << ": [" << desc.getX(i) << ", " << desc.getY(i) << ", " << desc.getZ(i) << "]" << "}  ";
  }
  return os;
}

}
