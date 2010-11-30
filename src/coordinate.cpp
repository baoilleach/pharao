/*******************************************************************************
coordinate.cpp - Pharao
 
Copyright (C) 2005-2010 by Silicos NV
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*******************************************************************************/



#include "coordinate.h"


		
Coordinate::Coordinate(void):
   x(0.0),
   y(0.0),
   z(0.0)
{
}


         
Coordinate::Coordinate(double x, double y, double z):
   x(x),
   y(y),
   z(z)
{
};



std::ostream&
operator<< (std::ostream& os, const Coordinate& A)
{
	os << "(" << A.x << "," << A.y << "," << A.z << ")";
	return os;
};
