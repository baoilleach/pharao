/*******************************************************************************
options.cpp - Pharao
 
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



#include "options.h"



Options::Options(): 
   cutOff(0.0),
   best(0),
   rankby(TANIMOTO),
   funcGroupVec(),
	noHybrid(false),
   epsilon(0.5),
   withExclusion(false), 
   merge(false),
   noNormal(false),
   scoreOnly(false),
   isQuiet(false),
   version(false),
   help(false)
{
   molOutStream = NULL;
   molOutWriter = NULL;
   molOutFile = "";

   pharmOutStream = NULL;
   pharmOutWriter = NULL;
   pharmOutFile = "";

   scoreOutStream = NULL;
   scoreOutFile = "";

   dbInpStream = NULL;
   dbInpFile = "";
	dbInpType = UNKNOWN;

   refInpStream = NULL;
   refInpFile = "";
   refInpType = UNKNOWN;
}




Options::~Options(void)
{
   //
   // reference input
   if (!refInpFile.empty())
   {
      refInpFile = "";
   };
   if (refInpStream)
   {
      delete refInpStream;
      refInpStream = NULL;
   };


   //
   // Database input
   if (!dbInpFile.empty())
   {
      dbInpFile = "";
   };
   if (dbInpStream)
   {
      delete dbInpStream;
      dbInpStream = NULL;
   };


   //
   // Molecule output
   if (!molOutFile.empty())
   {
      molOutFile = "";
   };
   if (molOutStream)
   {
      delete molOutStream;
      molOutStream = NULL;
   };
   if (molOutWriter)
   {
      delete molOutWriter;
      molOutWriter = NULL;
   };
   
   
   //
   // Pharmacopore output
   if (!pharmOutFile.empty())
   {
      pharmOutFile = "";
   };
   if (pharmOutStream)
   {
      delete pharmOutStream;
      pharmOutStream = NULL;
   };
   if (pharmOutWriter)
   {
      delete pharmOutWriter;
      pharmOutWriter = NULL;
   };
   
   
   //
   // Score output
   if (!scoreOutFile.empty())
   {
      scoreOutFile = "";
   };
   if (scoreOutStream)
   {
      delete scoreOutStream;
      scoreOutStream = NULL;
   };
}






std::string
Options::print(void) const
{
   std::ostringstream os;
   os << std::endl;
   os << "COMMAND_LINE OPTIONS:" << std::endl;
   os << std::endl;
   os << "  -> Reference file:    " << (refInpFile.empty() ? "no" : refInpFile) << std::endl;
   os << "  -> Reference type:    ";
   if (refInpType == MOL)
   {
      os << "MOL" << std::endl;
   }
   else if (refInpType == PHAR)
   {
      os << "PHAR" << std::endl;
   }
   else
   {
      os << "UNKNOWN" << std::endl;
   }
   os << "  -> Database file:     " << (dbInpFile.empty() ? "no" : dbInpFile) << std::endl;
   os << "  -> Database type:     ";
   if (dbInpType == MOL)
   {
      os << "MOL" << std::endl;
   }
   else if (dbInpType == PHAR)
   {
      os << "PHAR" << std::endl;
   }
   else
   {
      os << "UNKNOWN" << std::endl;
   }
   os << "  -> Mol output file:   " << (molOutFile.empty() ? "no" : molOutFile) << std::endl;
   os << "  -> Pharm output file: " << (pharmOutFile.empty() ? "no" : pharmOutFile) << std::endl;
   os << "  -> Scores file:       " << (scoreOutFile.empty() ? "no" : scoreOutFile) << std::endl;
   os << "  -> Cutoff:            ";
   if (cutOff)
   {
      os << cutOff << std::endl;
   }
   else
   {
      os << "no" << std::endl;
   }
   os << "  -> Best hits:         ";
   if (best)
   {
      os << best << std::endl;
   }
   else
   {
      os << "no" << std::endl;
   }
   os << "  -> Rank by:           " << rankby << std::endl;
   os << "  -> Functional groups: ";
   if (funcGroupVec[AROM]) os << "AROM ";
   if (funcGroupVec[HDON]) os << "HDON ";
   if (funcGroupVec[HACC]) os << "HACC ";
   if (funcGroupVec[LIPO]) os << "LIPO ";
   if (funcGroupVec[NEGC]) os << "NEGC ";
   if (funcGroupVec[POSC]) os << "POSC ";
   if (funcGroupVec[HYBH]) os << "HYBH ";
   if (funcGroupVec[HYBL]) os << "HYBL ";
   os << std::endl;
   os << "  -> Hybrids:           " << (noHybrid ? "no" : "yes") << std::endl;
   os << "  -> Epsilon:           " << epsilon << std::endl;
   os << "  -> Merge pharm:       " << (merge ? "yes" : "no") << std::endl;
   os << "  -> Include normals:   " << (noNormal ? "no" : "yes") << std::endl;
   os << "  -> With exclusion:    " << (withExclusion ? "yes" : "no") << std::endl;
   os << "  -> Scores only:       " << (scoreOnly ? "yes" : "no") << std::endl;
   os << "  -> Quied mode:        " << (isQuiet ? "yes" : "no") << std::endl;
   
   os << std::endl;   
   std::string r = os.str();
   return r;
}
