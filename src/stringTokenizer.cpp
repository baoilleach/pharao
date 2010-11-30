/*******************************************************************************
stringTokenizer.cpp - Pharao
 
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



#include "stringTokenizer.h"



std::list<std::string>
stringTokenizer(const std::string& s, const std::string& spacer)
{
   const std::string::size_type len = s.length();
   std::string::size_type i = 0;
   std::list<std::string> container;
   container.clear();
			
   while (i < len)
   {
      // eat leading whitespace
      i = s.find_first_not_of(spacer, i);
      if (i == std::string::npos)
      {
         return container;   // nothing left but white space
      }
				
      // find the end of the token
      std::string::size_type j = s.find_first_of(spacer, i);
				
      // push token onto container
      if (j == std::string::npos) // end of string
      {
         container.push_back(s.substr(i));
         return container;
      } 
      else // 
      {
         container.push_back(s.substr(i, j-i));
      }
      // move positions to next substring
      i = j + 1;
   }
   
   return container;
};
