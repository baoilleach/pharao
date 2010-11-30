/*******************************************************************************
utilities.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_UTILITIES_H__
#define __SILICOS_PHARAO_UTILITIES_H__



// General

// OpenBabel
#include "openbabel/mol.h"

// Pharao
#include "coordinate.h"
#include "pharmacophore.h"
#include "siMath.h"
#include "solutionInfo.h"

#ifndef GCI
#define GCI 2.828427125
#endif

#ifndef GCI2
#define GCI2 7.999999999
#endif

#ifndef PI
#define PI 3.14159265
#endif



Coordinate translate(Coordinate& p, Coordinate& t);
Coordinate rotate(Coordinate& p, SiMath::Matrix& R);
void normalise(Coordinate & p);
double norm(Coordinate & p);
double dotProduct(Coordinate& p1, Coordinate& p2) ;
Coordinate crossProduct(Coordinate& p1, Coordinate& p2);
double cosine(Coordinate& p1, Coordinate& p2);
double distance(Coordinate& p1, Coordinate& p2);
		
void normalise(SiMath::Vector& v);
		
SiMath::Matrix quat2Rotation(SiMath::Vector& Q);
		
void inverseHessian(SiMath::Matrix& H);

double VolumeOverlap(PharmacophorePoint& p1, PharmacophorePoint& p2, bool n);
double VolumeOverlap(PharmacophorePoint* p1, PharmacophorePoint* p2, bool n);
		
void TransformPharmacophore(Pharmacophore& pharm, SiMath::Matrix& U, Coordinate& center1, Coordinate& center2);
void positionPharmacophore(Pharmacophore& pharm, SiMath::Matrix& U, SolutionInfo& s);

void TransformMolecule(OpenBabel::OBMol* m, SiMath::Matrix& U, Coordinate& center1, Coordinate& center2);
void positionMolecule(OpenBabel::OBMol* m, SiMath::Matrix& U, SolutionInfo& s);
		


#endif

