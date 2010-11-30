/*******************************************************************************
hDonFuncCalc.cpp - Pharao
 
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



#include "hDonFuncCalc.h"



void
hDonFuncCalc(OpenBabel::OBMol* mol, Pharmacophore* pharmacophore)
{
   // Create for every hydrogen donor a pharmacophore point
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = mol->BeginAtom(ai); a; a = mol->NextAtom(ai))
   {
      if (a->GetAtomicNum() == 7 || a->GetAtomicNum() == 8)
      {
         if (a->GetFormalCharge() >= 0 && ((a->GetImplicitValence() - a->GetHvyValence()) !=0)) 
         {
            PharmacophorePoint p;
            p.func = HDON;
            p.point.x = a->x();
            p.point.y = a->y();
            p.point.z = a->z();
            p.hasNormal = true;
            p.alpha = funcSigma[HDON];
            p.normal = _hDonCalcNormal(a);
            pharmacophore->push_back(p);
         }
      }
   }
}



Coordinate
_hDonCalcNormal(OpenBabel::OBAtom* a)
{
   int nbrBonds(0);
   Coordinate normal;

   std::vector<OpenBabel::OBBond*>::iterator bi;
   for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
   {
      OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
      if (aa->GetAtomicNum() == 1)
      {
         continue;
      }
      ++nbrBonds;
      normal.x += (aa->x() - a->x());
      normal.y += (aa->y() - a->y());
      normal.z += (aa->z() - a->z());
   }
   double length(sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z));
   normal.x /= length;
   normal.y /= length;
   normal.z /= length;
  
   normal.x = -normal.x;
   normal.y = -normal.y;
   normal.z = -normal.z;
  
   normal.x += a->x();
   normal.y += a->y();
   normal.z += a->z();
  
   return normal;
}
