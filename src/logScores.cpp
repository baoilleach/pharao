/*******************************************************************************
logScores.cpp - Pharao
 
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



#include "logScores.h"



void
logScores(Result* res, Options& uo)
{
   *uo.scoreOutStream  << res->refId << "\t" 
	<< std::setprecision(6) << res->refVolume << "\t" 
	<< res->dbId << "\t" 
	<< std::setprecision(6) <<  res->dbVolume  << "\t" 
	<< std::setprecision(6) <<  res->overlapVolume << "\t"
	<< std::setprecision(6) <<  res->exclVolume << "\t"
	<< std::setprecision(6) <<  (res->info).volume << "\t" 
	<< res->resPharSize <<  "\t" 
	<< std::setprecision(4) <<  res->tanimoto << "\t"
	<< std::setprecision(4) <<  res->tversky_ref << "\t"
	<< std::setprecision(4) <<  res->tversky_db << std::endl;
}
