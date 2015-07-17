//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   In vrs 2.3.1, most of the driver implementation code was factored out 
   to the new GAtmoFlux base class so as to be shared by the functionally 
   identical GBartolAtmoFlux driver
 @ Feb 05, 2008 - Chris Backhouse (Oxford)
   Fixed a bug in bin definitions (the TH2 constructor takes an array of 
   bin lower edges, but the bin centres were being passed in instead).
 @ Feb 23, 2010 - CA
   Build bin arrays at ctor. Re-structuring and clean-up.

*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>
#include <string>

#include <TH3D.h>
#include <TMath.h>
#include <TFile.h>


#include "FluxDrivers/GHondaAtmoFlux.h"
#include "Messenger/Messenger.h"

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GHondaAtmoFlux::GHondaAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
       << "Instantiating the Honda atmospheric neutrino flux driver";

  this->SetBinSizes();
  this->Initialize();
}
//___________________________________________________________________________
GHondaAtmoFlux::~GHondaAtmoFlux()
{

}
//___________________________________________________________________________
void GHondaAtmoFlux::SetBinSizes(void)
{
// Generate the correct cos(theta), energy and phi bin sizes
// The flux is given in 20 bins of cos(zenith angle) from -1.0 to 1.0
// (bin width = 0.1), 12 bins of phi angle from 0 to 2Pi (bin width = Pi/6) 
// and 61 equally log-spaced energy bins (20 bins 
// per decade), with Emin = 0.100 GgeV.
//

  fCosThetaBins  = new double [kGHondaNumCosThetaBins  + 1];
  fEnergyBins    = new double [kGHondaNumLogEvBins     + 1];
  fPhiBins       = new double [kGHondaNumPhiBins       + 1];

  double dcostheta = 
      (kGHondaCosThetaMax - kGHondaCosThetaMin) /
      (double) kGHondaNumCosThetaBins;

  double dphi = 
    (kGHondaPhiMax - kGHondaPhiMin) /
    (double) kGHondaNumPhiBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(kGHondaEvMin);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGHondaNumLogEvBinsPerDecade;

  for(unsigned int i=0; i<= kGHondaNumCosThetaBins; i++) {
     fCosThetaBins[i] = kGHondaCosThetaMin + i * dcostheta;
     if(i != kGHondaNumCosThetaBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << kGHondaNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGHondaNumCosThetaBins];
     }
  }

  for(unsigned int i=0; i<= kGHondaNumPhiBins; i++) {
     fPhiBins[i] = kGHondaPhiMin + i * dphi;
     if(i != kGHondaNumPhiBins) {
       LOG("Flux", pDEBUG) 
         << "Honda flux: Phi bin " << i+1 
         << ": lower edge = " << fPhiBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "Honda flux: Phi bin " << kGHondaNumPhiBins 
         << ": upper edge = " << fPhiBins[kGHondaNumPhiBins];
     }
  }

  for(unsigned int i=0; i<= kGHondaNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kGHondaNumLogEvBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1 
         << ": lower edge = " << fEnergyBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << kGHondaNumLogEvBins 
         << ": upper edge = " << fEnergyBins[kGHondaNumLogEvBins];
     }
  }

  for(unsigned int i=0; i< kGHondaNumLogEvBins; i++) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": bin centre = " << (fEnergyBins[i] + fEnergyBins[i+1])/2.;
  }

  fNumCosThetaBins = kGHondaNumCosThetaBins;
  fNumEnergyBins   = kGHondaNumLogEvBins;
  fNumPhiBins      = kGHondaNumPhiBins;
}
//____________________________________________________________________________
bool GHondaAtmoFlux::FillFluxHisto3D(TH3D * histo, string filename, const int& pdg_nu)
{
  LOG("Flux", pNOTICE) << "Loading: " << filename;

  if(!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
  }
  ifstream flux_stream(filename.c_str(), ios::in);
  if(!flux_stream) {
     LOG("Flux", pERROR) << "Could not open file: " << filename;
     return false;
  }

  int    ibin, section, subsection, line;
  double energy, costheta, flux, phi;
  std::string junk;
  section = subsection = line = 1; //initialising some values
  costheta= 0.95;
  phi = 2.0*TMath::Pi()*(15.0/360.0); 

  double scale = 1.0; // 1.0 [m^2], OR 1.0e-4 [cm^2]

  if(pdg_nu == 14){
    while ( flux_stream ) {
      flux = 0.0;
      if (line == 1 || line == 2){
        std::getline(flux_stream, junk);
        line++; //ignore these lines
      } else {
        flux_stream >> energy >> flux >> junk >> junk >> junk; 
        line++;
        costheta = 1.0 -((double)section*0.1) + 0.05; 
        phi = 2.0*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0); 
        if( line == 104 ){ //new phi range
          ++subsection;
          line = 1;
          getline(flux_stream, junk);
          if (subsection == 13) //new costheta range
          {
            ++section;
            subsection = 1;
          }
        }
      }
      if( flux>0.0 ){
        LOG("Flux", pINFO)
          << "Flux[Ev = " << energy 
          << ", cos8 = " << costheta
          << ", phi = " << phi << "] = " << flux;
        ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta), (Axis_t)phi );   
        histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
      }
    }
  } else if(pdg_nu == -14){
    while ( flux_stream ) {
      flux = 0.0;
      if (line == 1 || line == 2){
        std::getline(flux_stream, junk);
        line++; //ignore these lines
      } else {
        flux_stream >> energy >> junk >> flux >> junk >> junk;
        line++;
        costheta = 1.0 -((double)section*0.1) + 0.05; 
        phi = 2*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0);
        if( line == 104 ){ //new phi range
          ++subsection;
          line = 1;
          getline(flux_stream, junk);
          if (subsection == 13) //new costheta range
          {
            ++section;
            subsection = 1;
          }
        }
      }
      if( flux>0.0 ){
        LOG("Flux", pINFO)
          << "Flux[Ev = " << energy 
          << ", cos8 = " << costheta
          << ", phi = " << phi << "] = " << flux;
        ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta), (Axis_t)phi );   
        histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
      }
    }
  } else if(pdg_nu == 12){
    while ( flux_stream ) {
      flux = 0.0;
      if (line == 1 || line == 2){
        std::getline(flux_stream, junk);
        line++; //ignore these lines
      } else {
        flux_stream >> energy >> junk >> junk >> flux >> junk; 
        line++;
        costheta = 1.0 -((double)section*0.1) + 0.05; 
        phi = 2*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0); 
        if( line == 104 ){ //new phi range
          ++subsection;
          line = 1;
          getline(flux_stream, junk);
          if (subsection == 13) //new costheta range
          {
            ++section;
            subsection = 1;
          }
        }
      }
      if( flux>0.0 ){
        LOG("Flux", pINFO)
          << "Flux[Ev = " << energy 
          << ", cos8 = " << costheta
          << ", phi = " << phi << "] = " << flux;
        ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta), (Axis_t)phi );   
        histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
      }
    }
  } else if (pdg_nu == -12){
    while ( flux_stream ) {
      flux = 0.0;
      if (line == 1 || line == 2){
        std::getline(flux_stream, junk);
        line++; //ignore these lines
      } else {
        flux_stream >> energy >> junk >> junk >> junk >> flux; 
        costheta = 1.0 -((double)section*0.1) + 0.05; 
        phi = 2*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0);  
        if( line == 104 ){ //new phi range
          ++subsection;
          line = 1;
          getline(flux_stream, junk);
          if (subsection == 13){ //new costheta range
            ++section;
            subsection = 1;
          }
        }
      }
      if( flux>0.0 ){
        LOG("Flux", pINFO)
          << "Flux[Ev = " << energy 
          << ", cos8 = " << costheta
          << ", phi = " << phi << "] = " << flux;
        ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta), (Axis_t)phi );   
        histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
      }
    }
  } else {
    LOG("FLUX", pERROR) 
      << "PDG code is not a neutrino type supported by this file.";
  }

  TFile f("starthisto.root","new");
  histo->Write();

  return true;
}
//___________________________________________________________________________
