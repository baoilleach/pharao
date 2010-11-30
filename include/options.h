/*******************************************************************************
options.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_OPTIONS_H__
#define __SILICOS_PHARAO_OPTIONS_H__



// General
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// OpenBabel
#include "openbabel/obconversion.h"

// Pharao
#include "fileType.h"
#include "rankType.h"
#include "pharmacophore.h"



class Options
{
   public:
   
      std::string             refInpFile;       //  -r  --reference
      FileType                refInpType;       //      --refType
      std::ifstream*          refInpStream;
      
      std::string             dbInpFile;        //  -d  --dbase
      FileType                dbInpType;        //      --dbType
      std::ifstream*          dbInpStream;
	
      std::string             pharmOutFile;     //  -p  --pharmacophore
      std::ofstream*          pharmOutStream;
      PharmacophoreWriter*    pharmOutWriter;
      
      std::string             molOutFile;       //  -o  --out
      std::ofstream*          molOutStream;
      OpenBabel::OBConversion* molOutWriter;
      
      std::string             scoreOutFile;     //  -s  --scores
      std::ofstream*          scoreOutStream;
      
      double                  cutOff;           //      --cutOff
      int                     best;             //      --best
      RankType                rankby;           //      --rankby
  
      std::vector<bool>       funcGroupVec;     //  -f  --funcGroup
      bool                    noHybrid;         //      --noHybrid
      double                  epsilon;          //  -e  --epsilon
      bool							withExclusion;    //      --withExclusion
      bool                    merge;            //  -m  --merge
      bool                    noNormal;         //  -n  --noNormal
      bool                    scoreOnly;        //      --scoreOnly
  
      bool                    isQuiet;          //  -q  --quiet
      bool                    version;          //  -v  --version
      bool                    help;             //  -h  --help
   
      Options(void);
      ~Options(void);
      
      std::string print(void) const;
};



#endif
