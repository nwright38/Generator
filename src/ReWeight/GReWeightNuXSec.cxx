//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "BaryonResonance/BaryonResonance.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNuXSec.h"
#include "ReWeight/GSystSet.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSec::GReWeightNuXSec() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSec::~GReWeightNuXSec()
{

}
//_______________________________________________________________________________________
void GReWeightNuXSec::Init(void)
{
  // Get the default cross section parameters 
  fXSecRwParams.LoadDefaults();
}
//_______________________________________________________________________________________
bool GReWeightNuXSec::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kSystNuXSec_MaQEL       ) : 
     case ( kSystNuXSec_MvQEL       ) : 
     case ( kSystNuXSec_MaRES       ) : 
     case ( kSystNuXSec_MvRES       ) : 
     case ( kSystNuXSec_MaCOHPi     ) : 
     case ( kSystNuXSec_RvpCC1pi    ) : 
     case ( kSystNuXSec_RvpCC2pi    ) : 
     case ( kSystNuXSec_RvpNC1pi    ) : 
     case ( kSystNuXSec_RvpNC2pi    ) : 
     case ( kSystNuXSec_RvnCC1pi    ) : 
     case ( kSystNuXSec_RvnCC2pi    ) : 
     case ( kSystNuXSec_RvnNC1pi    ) : 
     case ( kSystNuXSec_RvnNC2pi    ) : 
     case ( kSystNuXSec_RvbarpCC1pi ) : 
     case ( kSystNuXSec_RvbarpCC2pi ) :
     case ( kSystNuXSec_RvbarpNC1pi ) :
     case ( kSystNuXSec_RvbarpNC2pi ) : 
     case ( kSystNuXSec_RvbarnCC1pi ) : 
     case ( kSystNuXSec_RvbarnCC2pi ) : 
     case ( kSystNuXSec_RvbarnNC1pi ) : 
     case ( kSystNuXSec_RvbarnNC2pi ) : 

          handle = true;
          break;

     case ( kSystINuke_MFPTwk_pi    ) : 
     case ( kSystINuke_MFPTwk_N     ) : 
     case ( kSystINuke_CExTwk_pi    ) : 
     case ( kSystINuke_ElTwk_pi     ) : 
     case ( kSystINuke_InelTwk_pi   ) : 
     case ( kSystINuke_AbsTwk_pi    ) : 
     case ( kSystINuke_PiProdTwk_pi ) : 
     case ( kSystINuke_CExTwk_N     ) : 
     case ( kSystINuke_ElTwk_N      ) : 
     case ( kSystINuke_InelTwk_N    ) : 
     case ( kSystINuke_AbsTwk_N     ) : 
     case ( kSystINuke_PiProdTwk_N  ) : 
     default:

          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSec::SetSystematic(GSyst_t syst, double twk_dial)
{
   if( this->IsHandled(syst) ) {
      fXSecRwParams.SetCurTwkDial (syst, twk_dial);
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reset(void)
{
  fXSecRwParams.Reset();
  fXSecRwParams.Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reconfigure(void)
{
  fXSecRwParams.Reconfigure();
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcWeight(const genie::EventRecord & event) 
{
  if (! fXSecRwParams.IsTweaked() ) return 1.;

  double wght = fXSecRwHelper.NewWeight(event);
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcChisq()
{
  return fXSecRwParams.ChisqPenalty();
}
//_______________________________________________________________________________________