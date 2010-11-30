/*******************************************************************************
printHeader.cpp - Pharao
 
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



#include "printHeader.h"



void
printHeader()
{
	std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cerr << "  PHARAO v" << PHARAO_VERSION << "." << PHARAO_RELEASE << "." << PHARAO_SUBRELEASE << " | ";
   std::cerr << __DATE__ " " << __TIME__ << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -> GCC:         " << __VERSION__ << std::endl;
	std::cerr << "  -> Open Babel:  " << BABEL_VERSION << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Copyright (C) 2005-2010 by Silicos NV (http://www.silicos.com)" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  This program is part of the Open Babel project." << std::endl;
	std::cerr << "  For more information, see http://openbabel.sourceforge.net" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  This program is free software; you can redistribute it and/or modify" << std::endl;
	std::cerr << "  it under the terms of the GNU General Public License as published by" << std::endl;
	std::cerr << "  the Free Software Foundation version 2 of the License." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  This program is distributed in the hope that it will be useful," << std::endl;
	std::cerr << "  but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
	std::cerr << "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
	std::cerr << "  GNU General Public License for more details." << std::endl;
	std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cerr << std::endl;
}
