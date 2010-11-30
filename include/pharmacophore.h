/*******************************************************************************
pharmacophore.h - Pharao
 
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



#ifndef __SILICOS_PHARAO_PHARMACOPHORE_H__
#define __SILICOS_PHARAO_PHARMACOPHORE_H__



// General
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>

// OpenBabel

// Pharao
#include <coordinate.h>
#include <mainWar.h>
#include <siMath.h>
#include <stringTokenizer.h>

#ifndef GCI
#define GCI 2.828427125
#endif

#ifndef GCI2
#define GCI2 7.999999999
#endif

#ifndef PI
#define PI 3.14159265
#endif



enum FuncGroup
{
   AROM,  ///< Aromatic ringsystem, calculated by the AromFuncCalc class
   HDON,  ///< Hydrogen donor, calculated by the HDonFuncCalc class
   HACC,  ///< Hydrogen acceptor, calculated by the HAccFuncCalc class
   LIPO,  ///< Lipophilicity, calculated by the LipoFuncCalc class
   POSC,  ///< Positive charge, calculated by the ChargeFuncCalc class
   NEGC,  ///< Negative charge, calculated by the ChargeFuncCalc class
   HYBH,  ///< Hybrid Type: HDON + HACC
   HYBL,  ///< Hybrid Type: AROM + LIPO
   EXCL,  ///< Exclusion sphere
   UNDEF, ///< Undefined value (typically used for initialisation)
};



const std::string funcName[10] =
{
   "SILICOS::PHARAO::AROM",
   "SILICOS::PHARAO::HDON",
   "SILICOS::PHARAO::HACC",
   "SILICOS::PHARAO::LIPO",
   "SILICOS::PHARAO::POSC",
   "SILICOS::PHARAO::NEGC",
   "SILICOS::PHARAO::HYBH",
   "SILICOS::PHARAO::HYBL",
   "SILICOS::PHARAO::EXCL",
   "SILICOS::PHARAO::UNDEF"
};



const bool funcHasNormal[10] =
{
   true,
   true,
   true,
   false,
   false,
   false,
   true,
   true,
   false,
   false,
};



const double funcSigma[10] = 
{
   0.7, // AROM
   1.0, // HDON
   1.0, // HACC
   0.7, // LIPO
   1.0, // POSC
   1.0, // NEGC
   1.0, // HYBH
   0.7, // HYBL
   1.6, // EXCL
   1.0, // UNDEF
};



class PharmacophorePoint
{
   public:
   
      Coordinate     point;      ///< coordinates of the pharmacophore point
      FuncGroup      func;       ///< type of functional group
      Coordinate     normal;     ///< coordinates of the directionality of the pharmacophore
      double         alpha;      ///< spread of the gaussian is inverse proportional to radius squared
      bool           hasNormal;  ///< does the pharmacophore point has directionality
		
      PharmacophorePoint(void);
      PharmacophorePoint(const PharmacophorePoint&);
      PharmacophorePoint(const PharmacophorePoint*);
};



typedef std::vector<PharmacophorePoint> Pharmacophore;
typedef std::multimap<PharmacophorePoint*, PharmacophorePoint*> PharmacophoreMap;



class PharmacophoreReader
{
   private:

      void _skipPharmacophore(std::ifstream*);

   public:

      PharmacophoreReader(void);   
      ~PharmacophoreReader(void);
  
      Pharmacophore read(std::ifstream*, std::string&);
};



class PharmacophoreWriter
{
   public:
   
      PharmacophoreWriter(void);
      ~PharmacophoreWriter(void);
      
      void write(Pharmacophore&, std::ofstream*, const std::string&);
};
	


#endif