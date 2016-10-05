/*
* ArrayDescription.h
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
#ifndef __ARRAYDESCRIPTION_H
#define __ARRAYDESCRIPTION_H

#include <boost/tuple/tuple.hpp>

#include <map>

namespace mca {

      class ArrayDescription
      {
      public:

        typedef boost::tuple<double, double, double> ArrayPosition; /**< x, y, and z coordinates. */
        typedef int ElementId;

        ArrayDescription();
        virtual ~ArrayDescription();

        ElementId pushPosition(double x, double y, double z, const std::string &name = "");
        ElementId pushPosition(const ArrayPosition &position, const std::string &name = "");

        size_t size() const;
        bool empty() const;

        double distance(ElementId i, ElementId j) const;
        double distance(const std::string &i, const std::string &j) const;
        double maxDistance() const;
        double minDistance() const;

        void getPosition (ElementId id, ArrayPosition &position) const;
        void getPosition (const std::string &name, ArrayPosition &position) const;

        double getX(const ElementId &id) const;
        double getY(const ElementId &id) const;
        double getZ(const ElementId &id) const;

        void getX(std::vector<double> &x) const;
        void getY(std::vector<double> &y) const;
        void getZ(std::vector<double> &z) const;

        double getX(const std::string &name) const;
        double getY(const std::string &name) const;
        double getZ(const std::string &name) const;

        std::string getName(ElementId id) const;
        ElementId getId(const std::string &name) const;

        double getBandwidth() const;

        ArrayDescription static make_linear_array_description(const std::vector<double> &x);

      private:

        typedef std::map<std::string, int> ElementNames;
        typedef std::map<int, ArrayPosition> ElementList;
        typedef enum {X=0, Y=1, Z=2} AxisName;

        double get(const ElementId &id, AxisName axis) const;

        double distance(ElementList::const_iterator &i,
                        ElementList::const_iterator &j) const;

        ElementId getId(const ElementId &id) const;

        ElementNames _names;
        ElementList _list;

        ElementId getNewId();

      };

      std::ostream& operator<<(std::ostream& os, const ArrayDescription& desc);

}


#endif // ARRAYDESCRIPTION_H
