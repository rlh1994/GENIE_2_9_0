//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "CrossSections/COHXSecAR.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
COHXSecAR::COHXSecAR() :
XSecIntegratorI("genie::COHXSecAR")
{

}
//____________________________________________________________________________
COHXSecAR::COHXSecAR(string config) :
XSecIntegratorI("genie::COHXSecAR", config)
{

}
//____________________________________________________________________________
COHXSecAR::~COHXSecAR()
{

}
//____________________________________________________________________________
double COHXSecAR::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  const InitialState & init_state = in -> InitState();
  
  if(! model->ValidProcess(in) ) return 0.;
  
  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSecAR", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  
  Range1D_t y_lim = kps.Limits(kKVy);
  
  // Check this
  double Enu      = init_state.ProbeE(kRfLab);
  double Elep_min = (1.-y_lim.max) * Enu;
  double Elep_max = (1.-y_lim.min) * Enu;
  
  LOG("COHXSecAR", pINFO)
       << "Lepton energy integration range = [" << Elep_min << ", " << Elep_max << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);
  
  double xsec = 0;
  if (fSplitIntegral) {
    utils::gsl::dXSec_dElep_AR * func =
      new utils::gsl::dXSec_dElep_AR(model, interaction, fGSLIntgType, fGSLRelTol, fGSLMaxEval);
    
    //~ ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;
    ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
    
    double abstol = 1; // Pretty sure this parameter is unused by ROOT.
    int size = 1000;  // Max number of subintervals, won't reach nearly this.
    int rule = 2; // See https://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC283
                  // Rule 2 is 21 points min
    ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,size,rule);
    
    xsec = ig.Integral(Elep_min, Elep_max) * (1E-38 * units::cm2);
    delete func;
  }
  else {
    double zero    = kASmallNum;
    double pi      = kPi-kASmallNum ;
    double twopi   = 2*kPi-kASmallNum ;
    
    //~ ROOT::Math::IBaseFunctionMultiDim * func = 
          //~ new utils::gsl::wrap::d5Xsec_dEldOmegaldOmegapi(model, interaction);
    //~ double kine_min[5] = { Elep_min, zero , zero    , zero, zero };
    //~ double kine_max[5] = { Elep_max, pi   , twopi   , pi  , twopi};
    
    ROOT::Math::IBaseFunctionMultiDim * func = 
          new utils::gsl::d4Xsec_dEldThetaldOmegapi(model, interaction);
    double kine_min[4] = { Elep_min, zero , zero    , zero    };
    double kine_max[4] = { Elep_max, pi   , pi      , twopi   };
    
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
        
    double abstol = 1; //We mostly care about relative tolerance.
    ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  }

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void COHXSecAR::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::LoadConfig(void)
{
  // Get specified GENIE integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);

  // Get GSL integration type & relative tolerance
  fGSLIntgType   = fConfig->GetStringDef("gsl-integration-type" ,  "vegas");
  fGSLMaxEval    = (unsigned int) fConfig->GetIntDef("gsl-max-eval", 4000);
  fGSLRelTol     = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01); // Only used by adaptive
  fSplitIntegral = fConfig->GetBoolDef("split-integral", true);
}
//_____________________________________________________________________________


