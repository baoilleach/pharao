/*******************************************************************************
hAccFuncCalc.cpp - Pharao
 
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



#include "hAccFuncCalc.h"



void
hAccFuncCalc(OpenBabel::OBMol* mol, Pharmacophore* pharmacophore)
{
   // Create for every hydrogen acceptor a pharmacophore point
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* atom = mol->BeginAtom(ai); atom; atom = mol->NextAtom(ai))
   {
      if (atom->GetAtomicNum() == 7 || atom->GetAtomicNum() == 8)
      {
         if (atom->GetFormalCharge() <= 0)
         {        
            if(_hAccDelocalized(atom) || (_hAccCalcAccSurf(atom) < 0.02))
            {
               continue;
            }
            PharmacophorePoint p;
            p.func = HACC;
            p.point.x = atom->x();
            p.point.y = atom->y();
            p.point.z = atom->z();
            p.hasNormal = true;
            p.alpha = funcSigma[HACC];
            p.normal = _hAccCalcNormal(atom);
            pharmacophore->push_back(p);
         }
      }
   }
}



double
_hAccCalcAccSurf(OpenBabel::OBAtom* atom)
{
   double radius(H_BOND_DIST);
  
   //---(1)-- create sphere with uniformly distributed points
   std::vector<Coordinate> sphere;
   std::vector<Coordinate>::iterator itS;
  
   const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
   double dphi(arclength/radius);
	int nlayer(ROUND(PI/dphi) + 1);
  
   double phi(0.0);
   for (int i(0); i < nlayer; ++i) 
   {
      double rsinphi(radius*sin(phi));
      double z(radius*cos(phi));
      double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
      int tmpNbrPoints(ROUND(PI*2/dtheta));
      if(tmpNbrPoints <= 0)
      {
         tmpNbrPoints = 1;
      }
      dtheta = PI * 2.0 / tmpNbrPoints;
      double theta((i % 2) ? 0 : PI);
      for (int j(0) ; j < tmpNbrPoints ; ++j)
      {
         Coordinate coord;
         coord.x = rsinphi*cos(theta) + atom->x();
         coord.y = rsinphi*sin(theta) + atom->y();
         coord.z = z + atom->z();
         sphere.push_back(coord);
         theta += dtheta;
         if(theta > PI*2)
         {
            theta -= PI*2;
         }
      }
      phi += dphi;
   }

   //---(2)-- define neighbors of atom
   std::list<OpenBabel::OBAtom*> aList(_hAccGetNeighbors(atom));
   std::list<OpenBabel::OBAtom*>::iterator itA;
  
   //---(3) -- check for every sphere-point if it is accessible
   int nbrAccSurfPoints(0);
   double r;
   OpenBabel::OBElementTable et;
   for (itS = sphere.begin(); itS != sphere.end(); ++itS)
   {
      bool isAccessible(true);
      for (itA = aList.begin(); itA != aList.end(); ++itA) 
      {
         OpenBabel::OBAtom* n(*itA);
         double distSq(((itS->x - n->x()) * (itS->x - n->x())) +
                       ((itS->y - n->y()) * (itS->y - n->y())) +
                       ((itS->z - n->z()) * (itS->z - n->z())));
         radius = et.GetVdwRad(n->GetAtomicNum());
         double sumSq((r + H_RADIUS) * (r + H_RADIUS));
      
         if (distSq <= sumSq)
         {
            isAccessible = false;
            break;
         }
      }
    
      if (isAccessible)
      { 
         ++nbrAccSurfPoints;    
      
      }
   }
   
   return (nbrAccSurfPoints/(double)sphere.size());
}




std::list<OpenBabel::OBAtom*>
_hAccGetNeighbors(OpenBabel::OBAtom* a)
{
   std::list<OpenBabel::OBAtom*> aList;
   OpenBabel::OBMol* parent(a->GetParent());

   double r;
   OpenBabel::OBElementTable et;
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* aa = parent->BeginAtom(ai); aa; aa = parent->NextAtom(ai))
   {
      if (*aa == a)
      {
         continue;
      }
    
      r = et.GetVdwRad(aa->GetAtomicNum());
      double delta(H_BOND_DIST + H_RADIUS + r);
      double maxDistSq(delta*delta);
      double distSq((a->x() - aa->x()) * (a->x() - aa->x()) +
                    (a->y() - aa->y()) * (a->y() - aa->y()) +
                    (a->z() - aa->z()) * (a->z() - aa->z()));
    
      if (distSq <= maxDistSq)
      {
         aList.push_back(aa);
      }
   }
   return aList;
}



bool
_hAccDelocalized(OpenBabel::OBAtom* a)
{
   if (a->GetAtomicNum() != 7)
   {
      return false;
   }
   if (a->IsAromatic() && a->GetImplicitValence() == 3)
   {
      return true;
   }
  
   std::vector<OpenBabel::OBBond*>::iterator bi1;
   for (OpenBabel::OBBond* b1 = a->BeginBond(bi1); b1; b1 = a->NextBond(bi1))
   {
      OpenBabel::OBAtom* aa = b1->GetNbrAtom(a);
    
      if (aa->IsAromatic() && a->GetImplicitValence() == 3)
      {
         return true;
      }
    
      if (aa->GetAtomicNum() == 6)
      {
         std::vector<OpenBabel::OBBond*>::iterator bi2;
         for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2))
         {
            OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);
            
            if (aaa == a)
            {
               continue;
            }
            if (b2->GetBO() == 2)
            {
               if (aaa->GetAtomicNum() == 8)  return true;
               if (aaa->GetAtomicNum() == 7)  return true;
               if (aaa->GetAtomicNum() == 16) return true;
            }
         }
      }
      else if (aa->GetAtomicNum() == 16)
      {
         std::vector<OpenBabel::OBBond*>::iterator bi2;
         for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2))
         {
            OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);
            
            if (aaa == a)
            {
               continue;
            }
            if ((b2->GetBO() == 2) && (aaa->GetAtomicNum() == 8))
            {
               return true;
            }
         }
      }
   }
   return false;
}




Coordinate
_hAccCalcNormal(OpenBabel::OBAtom* a)
{
   Coordinate normal;
   std::vector<OpenBabel::OBBond*>::iterator bi;
   for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
   {
      OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
      if (aa->GetAtomicNum() == 1)
      {
         continue;
      }
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
