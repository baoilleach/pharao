/*******************************************************************************
parseCommandLine.cpp - Pharao
 
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



#include "parseCommandLine.h"



Options
parseCommandLine(int argc, char* argv[])
{
   static struct option Arguments[] =
   {
      { "reference",         required_argument,   NULL,   'r' },
		{ "dbase",             required_argument,   NULL,   'd' },
		{ "scores",            required_argument,   NULL,   's' },
		{ "out",               required_argument,   NULL,   'o' },
		{ "pharmacophore",     required_argument,   NULL,   'p' },
		{ "funcGroup",         required_argument,   NULL,   'f' },
		{ "epsilon",           required_argument,   NULL,   'e' },
      { "merge",             no_argument,         NULL,   'm' },
      { "noNormal",          no_argument,         NULL,   'n' },
		{ "help",              no_argument,         NULL,   'h' },
		{ "version",           no_argument,         NULL,   'v' },
		{ "quiet",             no_argument,         NULL,   'q' },
		{ "refType",           required_argument,   NULL,    1  },
		{ "dbType",            required_argument,   NULL,    2  },
      { "cutOff",            required_argument,   NULL,    3  },
      { "best",              required_argument,   NULL,    4  },
      { "rankBy",            required_argument,   NULL,    5  },
		{ "noHybrid",          no_argument,         NULL,    7  },
      { "info",              required_argument,   NULL,    9  },
		{ "withExclusion",     no_argument,         NULL,    10 },
		{ "scoreOnly",         no_argument,         NULL,    11 },
      { NULL,                0,                   NULL,    0  }
   };
	
	Options o;
	
	// Set defaults
	o.dbInpFile.clear();
   o.dbInpStream = NULL;
	o.dbInpType = UNKNOWN;

	o.refInpFile.clear();
   o.refInpStream = NULL;
	o.refInpType = UNKNOWN;
   
	o.molOutFile.clear();
   o.molOutStream = NULL;
   o.molOutWriter = NULL;
   
	o.pharmOutFile.clear();
   o.pharmOutStream = NULL;
   o.pharmOutWriter = NULL;
   
   o.scoreOutFile.clear();
   o.scoreOutStream = NULL;
   
	o.epsilon = 0.5;
  
   o.cutOff = 0.0;
   o.best = 0;
   o.rankby = TANIMOTO;
	
	o.funcGroupVec.resize(10);
	o.funcGroupVec[AROM] = true;
   o.funcGroupVec[HDON] = true;
	o.funcGroupVec[HACC] = true;
   o.funcGroupVec[LIPO] = true;
   o.funcGroupVec[POSC] = true;
   o.funcGroupVec[NEGC] = true;
   o.funcGroupVec[HYBH] = false;
   o.funcGroupVec[HYBL] = false;
	
	o.isQuiet = false;
   o.noHybrid = false;
   o.merge = false;
   o.noNormal = false;
	o.withExclusion = false;
	o.scoreOnly = false;
	o.version = false;
	
	int choice;
	opterr = 0;
	int optionIndex = 0;
	std::string strvalue, t;
	std::list<std::string> l;
	std::list<std::string>::iterator itL;
   std::string ext;
	
	while((choice = getopt_long(argc, argv,"vhqr:d:s:o:p:f:e:m", Arguments, &optionIndex )) != -1)
	{
      switch (choice)
		{
         case 'v': //....................................................version 
				o.version = true;
            break;
            
			case 'r': //..................................................reference 
            o.refInpFile = optarg;
            ext = getExt(o.refInpFile);
            if (ext == ".phar")
            {
               o.refInpType = PHAR;
            }
            else 
            {
               o.refInpType = MOL;
            }
            o.refInpStream = new std::ifstream(optarg);
            if (!o.refInpStream->good())
            {
               mainErr("Error opening input file for reference (-r)");
            }
            break;
            
			case 'd': //......................................................dbase
            o.dbInpFile = optarg;
            ext = getExt(o.dbInpFile);
            if (ext == ".phar")
            {
               o.dbInpType = PHAR;
            }
            else 
            {
               o.dbInpType = MOL;
            }
            o.dbInpStream = new std::ifstream(optarg);
            if (!o.dbInpStream->good())
            {
               mainErr("Error opening input file for database (-d)");
            }
            break;
            
			case 's': //.....................................................scores
            o.scoreOutFile = optarg;
            o.scoreOutStream = new std::ofstream(optarg);
            if (!o.scoreOutStream->good())
            {
               mainErr("Error opening output file for scores (-s)");
            }
            break;
            
			case 'o': //........................................................out
            o.molOutFile = optarg;
            o.molOutStream = new std::ofstream(optarg);
            if (!o.molOutStream->good())
            {
               mainErr("Error opening output file for molecules (-o)");
            }
            o.molOutWriter = new OpenBabel::OBConversion();
            o.molOutWriter->SetOutFormat(o.molOutWriter->FormatFromExt(optarg));
            break;
            
			case 'p': //..............................................pharmacophore
            o.pharmOutFile = optarg;
            o.pharmOutStream = new std::ofstream(optarg);
            if (!o.pharmOutStream->good())
            {
               mainErr("Error opening output file for pharmacophores (-p)");
            }
            o.pharmOutWriter = new PharmacophoreWriter();
            break;
            
			case 'e': //....................................................epsilon
            o.epsilon = strtod(optarg,NULL);
            break;
            
         case 'm': //......................................................merge
            o.merge = true;
            o.noNormal = true;
            break;
            
         case 'n': //...................................................noNormal
            o.noNormal = true;
            break;
            
			case 'f': //..................................................funcGroup
            {
               std::list<std::string> l = stringTokenizer(optarg, ",");
               std::list<std::string>::iterator itL;
               std::vector<bool> vec(10, false);
               for (itL = l.begin(); itL != l.end(); ++itL) 
               {
                  if (*itL == "AROM")
                  {
                     vec[AROM] = true;
                     continue;
                  }
                  if (*itL == "HDON")
                  {
                     vec[HDON] = true;
                     continue; 
                  }
                  if (*itL == "HACC")
                  {
                     vec[HACC] = true;
                     continue;
                  }
                  if (*itL == "LIPO")
                  {
                     vec[LIPO] = true;
                     continue;
                  }
                  if (*itL == "CHARGE")
                  {
                     vec[POSC] = true;
                     vec[NEGC] = true;
                     continue;
                  }
                  mainErr("Undefined functional Group. Only AROM, HDON, HACC, LIPO and"
                     "CHARGE are allowed as argument.");
               }
               o.funcGroupVec = vec; 
            }
            break;
         
			case 1: //......................................................refType
				t = optarg;
				if (t == "MOL")   
            {
               o.refInpType = MOL;
               break;
            }
            if (t == "PHAR")
            {
               o.refInpType = PHAR;
               break;
            }
				mainErr("Undefined reference type : " + t);
            break;
         
         case 2: //.......................................................dbType
				t = optarg;
				if (t == "MOL")
            {
               o.dbInpType = MOL;
               break;
            }
            if (t == "PHAR")
            {
               o.dbInpType = PHAR;
               break;
            }
            mainErr("Undefined dbase type: " + t);
				break;
            
         case 3: //.......................................................cutOff 
            o.cutOff = strtod(optarg, NULL);
            break;
         
         case 4: //.........................................................best 
            o.best = strtol(optarg, NULL, 10);
            break;
            
         case 5: //.......................................................rankby
            t = optarg;
            if (t == "TANIMOTO")
            {
               o.rankby = TANIMOTO;
               break;
            }
            if (t == "TVERSKY_REF")
            {
               o.rankby = TVERSKY_REF;
               break;
            }
            if (t == "TVERSKY_DB")
            {
               o.rankby = TVERSKY_DB;
               break;
            }
            mainErr("Undefined rankby type : " + t);
				break;
            
         case 7: //.....................................................noHybrid
				o.noHybrid = true;
            break;
            
         case 'h': //.......................................................help
            o.help = true;
            break;
            
         case 9: //.........................................................info
            printInfo(std::string(optarg));
            break;
         
         case 10: //...............................................withExclusion
            o.withExclusion = true;
            break;
				
         case 11: //...................................................scoreOnly
				o.scoreOnly = true;
            break;
            
         case 'q': //......................................................quiet
            o.isQuiet = true;
            break;
				
         default:
				mainErr("unknown command line option");
		}
	}
	
	// If no options are given print the help
	if (optind == 1) 
   {
		o.help = true;
	}
   
	argc -= optind;
	argv += optind;
	return o;
}
