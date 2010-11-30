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



#include "hybridCalc.h"



void
hybridCalc(OpenBabel::OBMol* m, Pharmacophore* pharmacophore)
{
   
   Pharmacophore::iterator itP;
   Pharmacophore::iterator itP2;
   for (itP = pharmacophore->begin(); itP != pharmacophore->end(); ++itP)
   {
      switch (itP->func)
      {
         //................................................................
         case AROM:
            itP2 = itP;
            ++itP2;
            for( ; itP2 != pharmacophore->end(); ++itP2) 
            {
               if (itP2->func == LIPO)
               {
                  if (_hybridSameHybLPoint(itP2->point, itP->point)) 
                  {
                     // modify first point
                     itP->func = HYBL;
                     itP->alpha = funcSigma[HYBL];
                     itP->point.x = (itP->point.x + itP2->point.x) / 2.0;
                     itP->point.y = (itP->point.y + itP2->point.y) / 2.0;
                     itP->point.z = (itP->point.z + itP2->point.z) / 2.0;
                     itP->hasNormal = false;
                     itP->normal.x = 0.0;
                     itP->normal.y = 0.0;
                     itP->normal.z = 0.0;
              
                     // remove second point
                     pharmacophore->erase(itP2);
                     --itP2;
                  }
               }
            }
            break;
      
         //................................................................
         case HACC:
            itP2 = itP;
            ++itP2;
            for( ; itP2 != pharmacophore->end(); ++itP2) 
            {
               if ((itP2->func == HDON) && (_hybridSameHybHPoint(itP2->point, itP->point)))
               {
                  //modify first point
                  itP->func = HYBH;
                  itP->alpha = funcSigma[HYBH];
                      
                  // Update normal
                  itP->normal.x = itP->normal.x - itP->point.x;
                  itP->normal.y = itP->normal.y - itP->point.y;
                  itP->normal.z = itP->normal.z - itP->point.z;
                     
                  itP2->normal.x = itP2->normal.x - itP2->point.x;
                  itP2->normal.y = itP2->normal.y - itP2->point.y;
                  itP2->normal.z = itP2->normal.z - itP2->point.z;
                     
                  itP->normal.x = (itP->normal.x + itP2->normal.x) / 2.0;
                  itP->normal.y = (itP->normal.y + itP2->normal.y) / 2.0;
                  itP->normal.z = (itP->normal.z + itP2->normal.z) / 2.0;
                     
                  double length(sqrt(itP->normal.x*itP->normal.x + 
                                     itP->normal.y*itP->normal.y + 
                                     itP->normal.z*itP->normal.z));
                     
                  itP->normal.x /= length;
                  itP->normal.y /= length;
                  itP->normal.z /= length;
                     
                  itP->normal.x += itP->point.x;
                  itP->normal.y += itP->point.y;
                  itP->normal.z += itP->point.z;
                                                    
                  //remove second point
                  pharmacophore->erase(itP2);
                  --itP2;
               }
            }
            break;

         //................................................................
         case HDON:
            itP2 = itP;
            ++itP2;
            for( ; itP2 != pharmacophore->end(); ++itP2) 
            {
               if (itP2->func == HACC)
               {
                  if (_hybridSameHybHPoint(itP2->point, itP->point)) 
                  {
                     // modify first point
                     itP->func = HYBH;
                     itP->alpha = funcSigma[HYBH];
                      
                     // Update normal
                     itP->normal.x = itP->normal.x - itP->point.x;
                     itP->normal.y = itP->normal.y - itP->point.y;
                     itP->normal.z = itP->normal.z - itP->point.z;
                     
                     itP2->normal.x = itP2->normal.x - itP2->point.x;
                     itP2->normal.y = itP2->normal.y - itP2->point.y;
                     itP2->normal.z = itP2->normal.z - itP2->point.z;
                     
                     itP->normal.x = (itP->normal.x + itP2->normal.x) / 2.0;
                     itP->normal.y = (itP->normal.y + itP2->normal.y) / 2.0;
                     itP->normal.z = (itP->normal.z + itP2->normal.z) / 2.0;
                     
                     double length(sqrt(itP->normal.x*itP->normal.x + 
                                        itP->normal.y*itP->normal.y + 
                                        itP->normal.z*itP->normal.z));
                     
                     itP->normal.x /= length;
                     itP->normal.y /= length;
                     itP->normal.z /= length;
                     
                     itP->normal.x += itP->point.x;
                     itP->normal.y += itP->point.y;
                     itP->normal.z += itP->point.z;
                                                    
                     //remove second point
                     pharmacophore->erase(itP2);
                     --itP2;
                  }
               }
            }
            break;

         //................................................................
         case LIPO:
            itP2 = itP;
            ++itP2;
            for( ; itP2 != pharmacophore->end(); ++itP2) 
            {
               if (itP2->func == AROM)
               {
                  if (_hybridSameHybLPoint(itP2->point, itP->point)) 
                  {
                     // Modify first point
                     itP->func = HYBL;
                     itP->alpha = funcSigma[HYBL];
                     itP->point.x = (itP->point.x + itP2->point.x) / 2.0;
                     itP->point.y = (itP->point.y + itP2->point.y) / 2.0;
                     itP->point.z = (itP->point.z + itP2->point.z) / 2.0;
                     itP->hasNormal = false;
                     itP->normal.x = 0.0;
                     itP->normal.y = 0.0;
                     itP->normal.z = 0.0;
              
                     // Remove second point
                     pharmacophore->erase(itP2);
                     --itP2;
                  }
               }
            }
            break;
      }
   }
      
   // For the lipophilic pharmacophores: rename both AROM and LIPO into HYBL
   for (itP = pharmacophore->begin(); itP != pharmacophore->end(); ++itP)
   {
      switch (itP->func)
      {
         //................................................................
         case AROM:
            itP->func = HYBL;
            itP->alpha = funcSigma[HYBL];
            itP->hasNormal = false;
            itP->normal.x = 0.0;
            itP->normal.y = 0.0;
            itP->normal.z = 0.0;
            break;
      
         //................................................................
         case LIPO:
            itP->func = HYBL;
            itP->alpha = funcSigma[HYBL];
            itP->hasNormal = false;
            itP->normal.x = 0.0;
            itP->normal.y = 0.0;
            itP->normal.z = 0.0;
            break;
      }
   }
}



bool
_hybridSameHybHPoint(const Coordinate& c1, const Coordinate& c2)
{
   double distSqr((c1.x - c2.x) * (c1.x - c2.x) +
                  (c1.y - c2.y) * (c1.y - c2.y) + 
                  (c1.z - c2.z) * (c1.z - c2.z));
   return distSqr < 0.0001;
}



bool
_hybridSameHybLPoint(const Coordinate& c1, const Coordinate& c2)
{
   double distSqr((c1.x - c2.x) * (c1.x - c2.x) +
                  (c1.y - c2.y) * (c1.y - c2.y) + 
                  (c1.z - c2.z) * (c1.z - c2.z));
   return distSqr < 1.0;
}

