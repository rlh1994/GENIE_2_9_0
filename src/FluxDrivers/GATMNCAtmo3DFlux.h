//____________________________________________________________________________
/*!

\class   genie::flux::GATMNCAtmo3DFlux

\brief   A driver for the 3-D atmospheric neutrino flux (commonly known as the
         `Honda flux') produced by the ATMNC cosmic ray simulation code.

\ref     Atmospheric neutrino flux at INO, South Pole and Pyhasalmi,
         M. Sajjad Athar (Aligarh Muslim U.) , M. Honda (Tokyo U., ICRR), 
         T. Kajita (Tokyo U., IPMU & Tokyo U., ICRR), K. Kasahara (Waseda U., RISE) 
         and S. Midorikawa (Aomori U.), Phys.Lett. B718 (2013) 1375-1380;
         http://inspirehep.net/record/1191405?ln=en

         The flux files necessary for running this flux driver can be obtained
         from:​http://www.icrr.u-tokyo.ac.jp/~mhonda/

\author  Ali Ajmi 
         Homi Bhabha National Institute, Mumbai

         Gobinda Majumder 
         Tata Institute of Fundamental Research, Mumbai

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 9, 2015

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         for the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_ATMNC_ATMO_3D_FLUX_I_H_
#define _GENIE_ATMNC_ATMO_3D_FLUX_I_H_

#include "FluxDrivers/GAtmoFlux.h"

namespace genie {
namespace flux  {

// Number of cos(zenith), azimuthat, and log(energy) bins in flux simulation

const unsigned int kGHnd3DNumCosThetaBins       = 40;
const double       kGHnd3DCosThetaMin           = -1.0;
const double       kGHnd3DCosThetaMax           = 1.0;
const unsigned int kGHnd3DNumPhiBins           	= 12;
const double       kGHnd3DPhiMin               	= 0.0;
const double       kGHnd3DPhiMax               	= 2.0*TMath::Pi();
const unsigned int kGHnd3DNumLogEvBins          = 101;	
const unsigned int kGHnd3DNumLogEvBinsPerDecade = 20;
const double       kGHnd3DEvMin                 = 0.095; // GeV

class GATMNCAtmo3DFlux: public GAtmoFlux {

public :
  GATMNCAtmo3DFlux();
 ~GATMNCAtmo3DFlux();

  //
  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files
  //

private:

  void SetBinSizes   (void);
  bool FillFluxHisto (TH3D * hist, string filename, const int& pdg_nu);
};

} // flux namespace
} // genie namespace

#endif // _GENIE_ATMNC_ATMO_3D_FLUX_I_H_
