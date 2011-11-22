/*******************************************************************************
lipoFuncCalc.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_LIPOFUNCCALC_H__
#define __SILICOS_PHARAO_LIPOFUNCCALC_H__



// General
#include <set>
#include <list>

// OpenBabel
#include "openbabel/mol.h"
#include "openbabel/atom.h"
#include "openbabel/bond.h"
#include "openbabel/ring.h"

// Pharao
#include "pharmacophore.h"





#ifndef ROUND
#define ROUND(x)      ((int) ((x) + 0.5))
#endif

#ifndef H_BOND_DIST
#define H_BOND_DIST   1.8
#endif

#ifndef H_RADIUS
#define H_RADIUS      1.2
#endif

#ifndef DENSITY
#define DENSITY       1.5
#endif

#ifndef PI
#define PI            3.14159265
#endif

#ifndef PROBE_RADIUS
#define PROBE_RADIUS  1.4
#endif

#ifndef REF_LIPO
#define REF_LIPO      9.87
#endif




void                          lipoFuncCalc(OpenBabel::OBMol*, Pharmacophore*);
void                          _lipoLabelAtoms(OpenBabel::OBMol*);
double                        _lipoCalcAccSurf(OpenBabel::OBAtom*);
void                          _lipoGroupAtoms(OpenBabel::OBMol*, Pharmacophore*);
std::list<OpenBabel::OBAtom*> _lipoGetNeighbors(OpenBabel::OBAtom*);
void                          _lipoLabelNeighbors(OpenBabel::OBAtom*, double);


#endif
