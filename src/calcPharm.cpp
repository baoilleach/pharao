/*******************************************************************************
calcPharm.cpp - Pharao
 
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



#include "calcPharm.h"



void
calcPharm(OpenBabel::OBMol* m, Pharmacophore* p, const Options& uo)
{
   if (uo.funcGroupVec[AROM]) { aromFuncCalc(m, p); }
   if (uo.funcGroupVec[HDON]) { hDonFuncCalc(m, p); }
   if (uo.funcGroupVec[HACC]) { hAccFuncCalc(m, p); }
   if (uo.funcGroupVec[LIPO]) { lipoFuncCalc(m, p); }
    
   if (uo.funcGroupVec[NEGC] || uo.funcGroupVec[POSC])
   {
      chargeFuncCalc(m, p); 
   }
   if (uo.funcGroupVec[HYBH] || uo.funcGroupVec[HYBL])
   {
      hybridCalc(m, p); 
   }
   return;
}
