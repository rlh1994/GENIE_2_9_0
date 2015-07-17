//____________________________________________________________________________
/*!

\class   genie::flux::GFlukaAtmo3DFlux

\brief   A flux driver for the FLUKA 3-D Atmospheric Neutrino Flux

\ref     Astrop.Phys.19 (2003) p.269; hep-ph/0207035; hep-ph/9907408
         Alfredo.Ferrari     <Alfredo.Ferrari@cern.ch>
         Paola.Sala          <Paola.Sala@cern.ch>
         Giuseppe Battistoni <Giuseppe.Battistoni@mi.infn.it>
         Teresa Montaruli    <Teresa.Montaruli@ba.infn.it>

         To be able to use this flux driver you will need to download the
         flux data from:  http://pcbat1.mi.infn.it/~battist/neutrino.html

         Please note that this class expects to read flux files formatted as 
         described in the above FLUKA flux page.
         Each file contains 3 columns:
	 - neutrino energy (GeV) at bin centre
	 - neutrino cos(zenith angle) at bin centre
         - neutrino flux (#neutrinos /GeV /m^2 /sec /sr)
         The flux is given in 40 bins of cos(zenith angle) from -1.0 to 1.0 
         (bin width = 0.05) and 61 equally log-spaced energy bins (20 bins per 
         decade), with Emin = 0.100 GeV.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created July 3, 2005 [during the most boring MINOS shift ever!]

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GHONDA_ATMO_FLUX_I_H_
#define _GHONDA_ATMO_FLUX_I_H_

#include "FluxDrivers/GAtmoFlux.h"
#include <TMath.h>

namespace genie {
namespace flux  {

// Number of cos(zenith) and energy bins in flux simulation
const unsigned int kGHondaNumCosThetaBins       = 20; 
const double       kGHondaCosThetaMin           = -1.0;
const double       kGHondaCosThetaMax           =  1.0;
const unsigned int kGHondaNumLogEvBins          = 101;
const unsigned int kGHondaNumLogEvBinsPerDecade = 20;
const double       kGHondaEvMin                 = 0.095; // GeV
const unsigned int kGHondaNumPhiBins            = 12; 
const double       kGHondaPhiMin                = 0.0;
const double       kGHondaPhiMax                = 2.0*TMath::Pi();

class GHondaAtmoFlux: public GAtmoFlux {

public :
  GHondaAtmoFlux();
 ~GHondaAtmoFlux();

  //
  // Most implementation is derived from the base GAtmoFlux
  // The concrete driver is only required to implement a function for
  // loading the input data files
  //

private:

  void SetBinSizes    (void);
  bool FillFluxHisto3D(TH3D * h3, string filename, const int& pdg_nu);
};

} // flux namespace
} // genie namespace

#endif // _GFLUKA_ATMO_3D_FLUX_I_H_
