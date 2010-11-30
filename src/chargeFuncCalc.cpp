/*******************************************************************************
chargeFuncCalc.cpp - Pharao
 
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



#include "chargeFuncCalc.h"



void
chargeFuncCalc(OpenBabel::OBMol* m, Pharmacophore* pharmacophore)
{
   // Create for every non-zero formal charge a pharmacophore point
   int charge;
   for(OpenBabel::OBMolAtomIter atom(m); atom; ++atom)
   {
      charge = atom->GetFormalCharge();
      if (charge < 0)
      {
         PharmacophorePoint p;
         p.func = NEGC;
         p.point.x = atom->x();
         p.point.y = atom->y();
         p.point.z = atom->z();
         p.alpha = funcSigma[NEGC];
         p.hasNormal = false;
         pharmacophore->push_back(p);
      }
      else if (charge > 0)
      {
         PharmacophorePoint p;
         p.func = POSC;
         p.point.x = atom->x();
         p.point.y = atom->y();
         p.point.z = atom->z();
         p.alpha = funcSigma[POSC];
         p.hasNormal = false;
         pharmacophore->push_back(p);
      }
   }
}