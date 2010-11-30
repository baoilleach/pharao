/*******************************************************************************
functionMapping.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_FUNCTIONMAPPING_H__
#define __SILICOS_PHARAO_FUNCTIONMAPPING_H__



// General
#include <vector>
#include <map>

// OpenBabel

// Pharao
#include "pharmacophore.h"
#include "utilities.h"



class FunctionMapping
{
   public:
      
      FunctionMapping(Pharmacophore*, Pharmacophore*, double);
      ~FunctionMapping(void);
		
      PharmacophoreMap getNextMap(void);
		
   private:
   
      bool _hasNext;
      unsigned int _maxLevel;
      double _epsilon;
			
      Pharmacophore* _ref;      ///< Pointer to reference pharmacophore
      Pharmacophore* _db;        ///< Pointer to database pharmacophore
						
      std::vector<unsigned int> _refIndex;
      std::vector<unsigned int> _dbIndex;
      std::map<unsigned int, std::vector<std::vector<unsigned int> > *> _matchMap;
};



#endif
