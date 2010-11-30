/*******************************************************************************
result.cpp - Pharao
 
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



#include "result.h"



Result::Result():
   refId(""),
   refVolume(0.0),
   dbId(""),
   dbVolume(0.0), 
   overlapVolume(0.0),
   exclVolume(0.0),
   resPharSize(0),
   tanimoto(0.0), 
   tversky_ref(0.0),
   tversky_db(0.0), 
   rankbyScore(0.0),
   info(),
   resMol(),
   resPhar()
{
};
