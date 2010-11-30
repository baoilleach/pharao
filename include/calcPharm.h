/*******************************************************************************
calcPharm.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_CALCPHARM_H__
#define __SILICOS_PHARAO_CALCPHARM_H__



// General

// OpenBabel
#include "openbabel/mol.h"

// Pharao
#include "pharmacophore.h"
#include "options.h"
#include "aromFuncCalc.h"
#include "chargeFuncCalc.h"
#include "hAccFuncCalc.h"
#include "hDonFuncCalc.h"
#include "lipoFuncCalc.h"
#include "hybridCalc.h"



void calcPharm(OpenBabel::OBMol*, Pharmacophore*, const Options&);




#endif
