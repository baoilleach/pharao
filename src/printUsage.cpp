/*******************************************************************************
printUsage.cpp - Pharao
 
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



#include "printUsage.h"



void
printUsage() 
{
   printHeader();
   std::cerr << "TASK: " << std::endl;
   std::cerr << std::endl;
   std::cerr << "  Pharao is a tool to generate and align pharmacophores." << std::endl;
   std::cerr << std::endl;
   std::cerr << "INPUT OPTIONS: " << std::endl;
   std::cerr << std::endl;
   std::cerr << "  -r  --reference       <file>" << std::endl;
   std::cerr << "      --refType         <MOL|PHAR>" << std::endl;
   std::cerr << "  -d  --dbase           <file>" << std::endl;
   std::cerr << "      --dbType          <MOL|PHAR>" << std::endl;
   std::cerr << std::endl;
   std::cerr << "OUTPUT OPTIONS: " << std::endl;
   std::cerr << std::endl;
   std::cerr << "  -p  --pharmacophore   <file>" << std::endl;
   std::cerr << "  -o  --out             <file>" << std::endl;
   std::cerr << "  -s  --scores          <file>" << std::endl;
   std::cerr << std::endl;
   std::cerr << "      --cutOff          <double>" << std::endl;
   std::cerr << "      --best            <int>" << std::endl;
   std::cerr << "      --rankBy          <TANIMOTO|TVERSKY_REF|TVERSKY_DB>" << std::endl;
   std::cerr << std::endl;  
   std::cerr << "PHARAO OPTIONS: " << std::endl;
   std::cerr << std::endl;
   std::cerr << "  -f  --funcGroup       <AROM|HDON|HACC|LIPO|CHARGE>" << std::endl;
   std::cerr << "  -e  --epsilon         <double>" << std::endl;
   std::cerr << "  -m  --merge" << std::endl;
   std::cerr << "  -n  --noNormal" << std::endl;
   std::cerr << "      --noHybrid" << std::endl;
   std::cerr << "      --withExclusion" << std::endl;
   std::cerr << "      --scoreOnly" << std::endl;
   std::cerr << std::endl;
   std::cerr << "GENERAL OPTIONS: " << std::endl;
   std::cerr << std::endl;
   std::cerr << "  -h  --help" << std::endl;
   std::cerr << "      --info            <option>" << std::endl;
   std::cerr << "  -q  --quiet" << std::endl;
   std::cerr << std::endl;
}
