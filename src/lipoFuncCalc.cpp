/*******************************************************************************
lipoFuncCalc.cpp - Pharao
 
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



#include "lipoFuncCalc.h"



void
lipoFuncCalc(OpenBabel::OBMol* m, Pharmacophore* pharmacophore)
{
   // Make a copy of the partial charge of each atom, since here we store the values
   std::vector<double> pcharge;
   pcharge.clear();
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai))
   {
      pcharge.push_back(a->GetPartialCharge());
   }

   // Find for each atom the 'topology-dependent term': t
   _lipoLabelAtoms(m);
  
   // Multiply 't' with the accessible surface 's'
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai))
   {
      if (a->GetAtomicNum() == 1)
      {
         continue;
      }
      double t = a->GetPartialCharge();
      if (t != 0.0)
      {
         a->SetPartialCharge(_lipoCalcAccSurf(a) * t);
      }
   }
  
   // Finally calculate the lipophilic points
   _lipoGroupAtoms(m, pharmacophore);
   
   // Reset partial charges
   unsigned int i(0);
   for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai))
   {
      atom->SetPartialCharge(pcharge[i]);
      ++i;
   }
}



void 
_lipoLabelAtoms(OpenBabel::OBMol* m)
{
  // Give all atoms default score of 1
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai))
   {
      atom->SetPartialCharge(1.0);
   }
  
   // Decrease scores based on O,N,S and charged atoms
   std::vector<OpenBabel::OBBond*>::iterator bi;
   std::vector<OpenBabel::OBBond*>::iterator bi2;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai))
   {
      switch (a->GetAtomicNum())
      {
         //.........................................................
         case 1: //.................................................
            a->SetPartialCharge(0.0); // category 1
            break;
        
         //.........................................................
         case 7: //.................................................
            a->SetPartialCharge(0.0); // category 1
            if (a->IsAromatic())
            {
               break;
            }
            _lipoLabelNeighbors(a, 0.25); // category 13
            if ((a->GetImplicitValence() - a->GetHvyValence()) !=0 )
            {
               std::vector<OpenBabel::OBBond*>::iterator bi;
               for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
               {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0);     // category 4
                  _lipoLabelNeighbors(aa, 0.0);  // category 4
               }
            }
            break;
        
         //.........................................................
         case 8: //.................................................
            a->SetPartialCharge(0.0); // category 1
            if (a->IsAromatic())
            {
               break;
            }
            _lipoLabelNeighbors(a, 0.25); // category 13
        
            for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
            {
               if (b->GetNbrAtom(a)->GetAtomicNum() == 1)
               {
                  std::vector<OpenBabel::OBBond*>::iterator bi2;
                  for (OpenBabel::OBBond* b2 = a->BeginBond(bi2); b2; b2 = a->NextBond(bi2))
                  {
                     OpenBabel::OBAtom* aa = b2->GetNbrAtom(a);
                     aa->SetPartialCharge(0.0);    // category 4
                     _lipoLabelNeighbors(aa, 0.0); // category 4
                  }
               }
               if (b->GetBO() == 2)
               {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0); // category 6
                  for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2))
                  {
                     OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);
                     if (aaa == a)
                     {
                        continue;
                     }
                     aaa->SetPartialCharge(0.0);    // category 6
                     _lipoLabelNeighbors(aaa, 0.6); // category 9
                  }
               }
            }
            break;
          
         //.........................................................
         case 16: //................................................
            for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
            {
               if (b->GetNbrAtom(a)->GetAtomicNum() == 1)
               {
                  a->SetPartialCharge(0.0);    // category 2
                  _lipoLabelNeighbors(a, 0.0); // category 5
               }
          
               if (b->GetBO() == 2) 
               {
                  a->SetPartialCharge(0.0);    // category 8
                  _lipoLabelNeighbors(a, 0.6); // category 11
               }
            }
          
            // we make the assumption that every S with valence > 2 can be found
            // by counting the number of bonds.
            if ((a->GetImplicitValence() - a->GetHvyValence()) > 2)
            {
               a->SetPartialCharge(0.0); // category 7
               for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
               {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0);    // category 7
                  _lipoLabelNeighbors(aa, 0.6); // category 10
               }
            }
            break;
      }
    
      if (a->GetFormalCharge() != 0)
      {
         for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
         {
            OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
            aa->SetPartialCharge(0.0);
            _lipoLabelNeighbors(aa, 0.0);
         }
      }
   }
  
   // Adjust combinations of scores (category 12 and 14)
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai))
   {
      double value = a->GetPartialCharge();
      if ((value == 0.36 || value < 0.25) && value != 0.15)
      {
         a->SetPartialCharge(0.0); // category 12 and 14
      }
   }
}
 


double
_lipoCalcAccSurf(OpenBabel::OBAtom* a)
{
   OpenBabel::OBElementTable et;
   double radius = et.GetVdwRad(a->GetAtomicNum());
  
   // Create sphere with uniformly distributed points
   std::vector<Coordinate> sphere;
   std::vector<Coordinate>::iterator itS;
  
   const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
   double dphi(arclength/radius);
	int nlayer(ROUND(PI/dphi) + 1);
  
   double phi(0.0);
   for (int i = 0; i < nlayer; ++i) 
   {
      double rsinphi(radius*sin(phi));
      double z(radius*cos(phi));
      double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
      int tmpNbrPoints(ROUND(PI*2/dtheta));
      if (tmpNbrPoints <= 0)
      {
         tmpNbrPoints = 1;
      }
      dtheta = PI*2 / tmpNbrPoints;
      double theta((i%2) ? 0 : PI);
      for (int j = 0; j < tmpNbrPoints ; ++j)
      {
         Coordinate coord;
         coord.x = rsinphi*cos(theta) + a->x();
         coord.y = rsinphi*sin(theta) + a->y();
         coord.z = z + a->z();
         sphere.push_back(coord);
         theta += dtheta;
         if (theta > PI*2)
         {
            theta -= PI*2;
         }
      }
      phi += dphi;
   }
  
   // Define neighbors of atom
   std::list<OpenBabel::OBAtom*> aList(_lipoGetNeighbors(a));
   std::list<OpenBabel::OBAtom*>::iterator itA;
  
  // Check for every sphere-point if it is accessible
  int nbrAccSurfPoints(0);
  double delta(PROBE_RADIUS/radius);
  for (itS = sphere.begin(); itS != sphere.end(); ++itS)
  {
      Coordinate p;
      p.x = ((itS->x - a->x()) * delta) + itS->x;
      p.y = ((itS->y - a->y()) * delta) + itS->y;
      p.z = ((itS->z - a->z()) * delta) + itS->z;
    
      bool isAccessible(true);
      for (itA = aList.begin(); itA != aList.end(); ++itA) 
      {
         OpenBabel::OBAtom* n(*itA);
         double distSq(((p.x - n->x()) * (p.x - n->x())) +
                       ((p.y - n->y()) * (p.y - n->y())) +
                       ((p.z - n->z()) * (p.z - n->z())));
         double sumSq((PROBE_RADIUS + et.GetVdwRad(n->GetAtomicNum())) *
                      (PROBE_RADIUS + et.GetVdwRad(n->GetAtomicNum())));

         if (distSq <= sumSq)
         { 
            isAccessible = false;
            break;
         }
      }
     
      if (isAccessible) ++nbrAccSurfPoints;    
   }
  
   double f(nbrAccSurfPoints/(double)sphere.size());
   return f * 4 * PI * radius*radius;
}
  


void 
_lipoGroupAtoms(OpenBabel::OBMol* m, Pharmacophore* pharmacophore)
{
   std::set<OpenBabel::OBAtom*> atomSet; // keeps remaining atoms for step (3)
   std::set<OpenBabel::OBAtom*>::iterator itS;
  
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai))
   {
      atomSet.insert(a);
   }
  
   // Group rings smaller than 7
   std::vector<OpenBabel::OBRing*> allrings = m->GetSSSR();
   std::vector<OpenBabel::OBRing*>::iterator ri;
   OpenBabel::OBRing* ring;
   for (ri = allrings.begin(); ri != allrings.end(); ++ri)
   {
      ring = *ri;
      if (ring->Size() > 7)
      {
         continue;
      }
      
      double lipoSum(0.0);
      Coordinate center;
      for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai))
      {
         if (ring->IsMember(atom))
         {
            atomSet.erase(atom);
            double lipo(atom->GetPartialCharge());
            lipoSum += lipo;
            center.x += lipo * atom->x();
            center.y += lipo * atom->y();
            center.z += lipo * atom->z();
         }
         else
         {
            continue;
         }
      }
         
      if (lipoSum > REF_LIPO)
      {
         PharmacophorePoint p;
         p.func = LIPO;
         p.hasNormal = false;
         p.normal.x = 0.0;
         p.normal.y = 0.0;
         p.normal.z = 0.0;
         p.alpha = funcSigma[LIPO];
         center.x /= lipoSum;
         center.y /= lipoSum;
         center.z /= lipoSum;
         p.point = center;
         pharmacophore->push_back(p);
      }
   }
  
   // Group atoms with three or more bonds
   for (itS = atomSet.begin(); itS != atomSet.end(); ++itS)
   {
      if (((*itS)->GetImplicitValence() - (*itS)->GetHvyValence()) > 2)
      {
         std::list<OpenBabel::OBAtom*> aList;
         aList.push_back(*itS);
         double lipoSum((*itS)->GetPartialCharge());
         Coordinate center;
         center.x += lipoSum * (*itS)->x();
         center.y += lipoSum * (*itS)->y();
         center.z += lipoSum * (*itS)->z();
         
         std::vector<OpenBabel::OBBond*>::iterator bi;
         for (OpenBabel::OBBond* b = (*itS)->BeginBond(bi); b; b = (*itS)->NextBond(bi))
         {
            OpenBabel::OBAtom* a = b->GetNbrAtom(*itS);
            if (((a->GetImplicitValence() - a->GetHvyValence()) == 1) && (a->GetAtomicNum() != 1))
            {
               double lipo(a->GetPartialCharge());
               lipoSum += lipo;
               aList.push_back(a);
               center.x += (lipo * a->x());
               center.y += (lipo * a->y());
               center.z += (lipo * a->z());
            }
         }
         
         if (lipoSum > REF_LIPO)
         {
            PharmacophorePoint p;
            p.func = LIPO;
            p.hasNormal = false;
            p.normal.x = 0.0;
            p.normal.y = 0.0;
            p.normal.z = 0.0;
				p.alpha = funcSigma[LIPO];
            center.x /= lipoSum;
            center.y /= lipoSum;
            center.z /= lipoSum;
            p.point = center;
            pharmacophore->push_back(p);
         }
      }
   }
}



std::list<OpenBabel::OBAtom*>
_lipoGetNeighbors(OpenBabel::OBAtom* a)
{
   OpenBabel::OBElementTable et;
   double radius = et.GetVdwRad(a->GetAtomicNum());
   std::list<OpenBabel::OBAtom*> aList;
  
   OpenBabel::OBMol* parent(a->GetParent());
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* aa = parent->BeginAtom(ai); aa; aa = parent->NextAtom(ai))
   {
      if ((aa->GetAtomicNum() == 1) || (aa == a))
      {
         continue;
      }
      
      double delta(radius + et.GetVdwRad(aa->GetAtomicNum()) + 2 * PROBE_RADIUS);
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



void
_lipoLabelNeighbors(OpenBabel::OBAtom* a, double value)
{ 
   std::vector<OpenBabel::OBBond*>::iterator bi;
   for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi))
   {
      OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
      aa->SetPartialCharge(value * aa->GetPartialCharge());
   }
}
