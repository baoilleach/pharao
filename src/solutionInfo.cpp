/*******************************************************************************
solutionInfo.cpp - Pharao
 
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



#include "solutionInfo.h"


		
SolutionInfo::SolutionInfo(void):
   rotor(4,0.0),
   volume(-999.99),
   iterations(0),
   center1(),
   center2(),
   rotation1(3,3,0.0),
   rotation2(3,3,0.0)
{
};


		
SolutionInfo::~SolutionInfo(void)
{
};
