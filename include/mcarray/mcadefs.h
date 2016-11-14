/*
* mcadefs.h
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
/*
* dspdefs.h
* Copyright 2016 (c) Jordi Adell
* Created on: 2015
* 	Author: Jordi Adell - adellj@gmail.com
*
* This file is part of DSPONE
*
* DSPONE is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPONE is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DSPONE.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ___MCADEFS_H_
#define ___MCADEFS_H_

#include <stdint.h>

#include <vector>
#include <mcarray/complex.h>

#include <boost/shared_array.hpp>

namespace mca{


// 32 bits float
typedef float BaseType32;
typedef Complex32f BaseType32C;

typedef boost::shared_array<BaseType32> SignalPtr32;
typedef boost::shared_array<BaseType32C> SignalCPtr32;
typedef std::vector<SignalPtr32> SignalVector32;
typedef std::vector<SignalCPtr32> SignalCVector32;


// 64 bits float
typedef double BaseType64;
typedef Complex64f BaseType64C;

typedef boost::shared_array<BaseType64> SignalPtr64;
typedef boost::shared_array<BaseType64C> SignalCPtr64;
typedef std::vector<SignalPtr64> SignalVector64;
typedef std::vector<SignalCPtr64> SignalCVector64;

// 16 bits integer
typedef signed short BaseType16s;
typedef boost::shared_array<BaseType16s> SignalPtr16s;
typedef std::vector<SignalPtr16s> SignalVector16s;


// Default types
typedef BaseType64 BaseType;
typedef BaseType64C BaseTypeC;
typedef boost::shared_array<BaseType> BaseTypePtr;
typedef boost::shared_array<BaseTypeC> BaseTypeCPtr;
typedef boost::shared_array<BaseType> SignalPtr;
typedef boost::shared_array<BaseTypeC> SignalCPtr;
typedef std::vector<SignalPtr> SignalVector;
typedef std::vector<SignalCPtr> SignalCVector;

}



#endif
