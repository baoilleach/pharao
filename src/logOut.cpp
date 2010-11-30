/*******************************************************************************
logOut.cpp - Pharao
 
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



#include "logOut.h"



void
logOut(Result* res, Options& uo)
{
   if((res->resMol).Empty())
   {
      return;
   }
   
   std::ostringstream ss;

	// Add properties to the molecule
   OpenBabel::OBPairData* label1 = new OpenBabel::OBPairData();
   label1->SetAttribute("PHARAO_TANIMOTO");
   ss << res->tanimoto;
   label1->SetValue(ss.str());
   (res->resMol).SetData(label1);
   
   OpenBabel::OBPairData* label2 = new OpenBabel::OBPairData();
   label2->SetAttribute("PHARAO_TVERSKY_REF");
   ss.str("");
   ss << res->tversky_ref;
   label2->SetValue(ss.str());
   (res->resMol).SetData(label2);
   
   OpenBabel::OBPairData* label3 = new OpenBabel::OBPairData();
   label3->SetAttribute("PHARAO_TVERSKY_DB");
   ss.str("");
   ss << res->tversky_db;
   label3->SetValue(ss.str());
   (res->resMol).SetData(label3);

   // Write molecule
   uo.molOutWriter->Write(&(res->resMol), uo.molOutStream);
   
   // Clean up
   (res->resMol).DeleteData("PHARAO_TANIMOTO");
   (res->resMol).DeleteData("PHARAO_TVERSKY_REF");
   (res->resMol).DeleteData("PHARAO_TVERSKY_DB");
}
