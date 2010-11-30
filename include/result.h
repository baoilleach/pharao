/*******************************************************************************
result.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_RESULT_H__
#define __SILICOS_PHARAO_RESULT_H__



// General
#include <string>
#include <vector>

// OpenBabel
#include <openbabel/mol.h>

// Pharao
#include "pharmacophore.h"
#include "solutionInfo.h"



class Result
{
   public:
   
      std::string             refId;         // id of the reference pharmacophore
      double                  refVolume;     // volume of the reference pharmacophore
      std::string             dbId;          // id of the database pharmacophore
      double                  dbVolume;      // volume of the database pharmacophore
      double                  overlapVolume; // volume overlap between reference and database pharmacophore
      double                  exclVolume;    // volume overlap between database pharmacophore and exclusion spheres
      int                     resPharSize;   // number of points in the resulting pharmacophore
	
      double                  tanimoto;      // resulting score = info.volume/(refVolume+resVolume-info.volume)
      double                  tversky_ref;   // info.volume/refVolume
      double                  tversky_db;    // info.volume/dbVolume
      double                  rankbyScore;   // one of the three scores...
  
      SolutionInfo            info;          // information about the alignment
      OpenBabel::OBMol        resMol;        // resulting molecule
      Pharmacophore           resPhar;       // overlapping pharmacophore
      
      Result(void);
};



#endif
