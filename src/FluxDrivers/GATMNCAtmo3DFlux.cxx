//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Ali Ajmi 
         Homi Bhabha National Institute, Mumbai

         Gobinda Majumder
         Tata Institute of Fundamental Research, Mumbai

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>
#include <string>

#include <TH3D.h>
#include <TMath.h>

#include "FluxDrivers/GATMNCAtmo3DFlux.h"
#include "Messenger/Messenger.h"

#include "FluxDrivers/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GATMNCAtmo3DFlux,genie::flux::GATMNCAtmo3DFlux)

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GATMNCAtmo3DFlux::GATMNCAtmo3DFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
    << "Instantiating the ATMNC 3D atmospheric neutrino flux driver";

  this->SetBinSizes();
  this->Initialize();
}
//___________________________________________________________________________
GATMNCAtmo3DFlux::~GATMNCAtmo3DFlux()
{

}
//___________________________________________________________________________
void GATMNCAtmo3DFlux::SetBinSizes(void)
{ 
  //
  // cos(theta)
  //

  fCosThetaBins    = new double [kGHnd3DNumCosThetaBins + 1];
  fNumCosThetaBins = kGHnd3DNumCosThetaBins;

  double dcostheta = 
      (kGHnd3DCosThetaMax - kGHnd3DCosThetaMin) /
      (double) kGHnd3DNumCosThetaBins;

  for(unsigned int i=0; i<= kGHnd3DNumCosThetaBins; i++) {
     fCosThetaBins[i] = kGHnd3DCosThetaMin + i * dcostheta;
     if(i != kGHnd3DNumCosThetaBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: CosTheta bin " << kGHnd3DNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGHnd3DNumCosThetaBins];
     }
  }


  //
  // phi
  //

  fPhiBins    = new double [kGHnd3DNumPhiBins + 1];
  fNumPhiBins = kGHnd3DNumPhiBins;

  double dphi = 
      (kGHnd3DPhiMax - kGHnd3DPhiMin) /
      (double) kGHnd3DNumPhiBins;

  for(unsigned int i=0; i<= kGHnd3DNumPhiBins; i++) {
     fPhiBins[i] = kGHnd3DPhiMin + i * dphi;
     if(i != kGHnd3DNumPhiBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Phi bin " << i+1 
         << ": lower edge = " << fPhiBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Phi bin " << kGHnd3DNumPhiBins 
         << ": upper edge = " << fPhiBins[kGHnd3DNumPhiBins];
     }
  }

  //
  // log(E)
  //

  fEnergyBins    = new double [kGHnd3DNumLogEvBins + 1];
  fNumEnergyBins = kGHnd3DNumLogEvBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(kGHnd3DEvMin);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGHnd3DNumLogEvBinsPerDecade;

  for(unsigned int i=0; i<= kGHnd3DNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kGHnd3DNumLogEvBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Energy bin " << i+1 
         << ": lower edge = " << fEnergyBins[i] << " GeV"
         << ", bin centre = " << (fEnergyBins[i] + fEnergyBins[i+1])/2. << " GeV";
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Energy bin " << kGHnd3DNumLogEvBins 
         << ": upper edge = " << fEnergyBins[kGHnd3DNumLogEvBins] << " GeV";
     }
  }

}
//____________________________________________________________________________
bool GATMNCAtmo3DFlux::FillFluxHisto(TH3D * histo, string filename)
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
        costheta = 1.0 -((double)section*0.1) + 0.05; //costheta is known based on what
                                                     //section of data we are in, this gives middle value
        phi = 2.0*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0); //phi known by subsection, again gives middle value
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
        line++;
        costheta = 1.0 -((double)section*0.1) + 0.05; //costheta is known based on what
                                                     //section of data we are in, this gives middle value
        phi = 2*TMath::Pi()*((-15.0 + ((double)subsection * 30.0))/360.0);  //phi known by subsection, again gives middle value
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
  return true;
}
//___________________________________________________________________________
