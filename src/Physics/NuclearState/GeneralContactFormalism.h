//____________________________________________________________________________
/*!

\class    genie::GeneralContactFormalism

\brief    local Fermi gas model. Implements the NuclearModelI 
          interface.

\ref      

\author   Joe Johnston, Steven Dytman

\created  December 2015

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _GENERALCONTACTFORMALISM_H_
#define _GENERALCONTACTFORMALISM_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"


namespace genie {

class GeneralContactFormalism : public NuclearModelI {

public:
  GeneralContactFormalism();
  GeneralContactFormalism(string config);
  virtual ~GeneralContactFormalism();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- allow methods to be called with a radius
  bool   GenerateNucleon (const Target & t, double hitNucleonRadius) const;
  double Prob            (double p, double w, const Target & t,
			  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmGeneralContactFormalism; 
  }

  virtual double LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set) ;

 protected:
  void   LoadConfig (void);
  TH1D * ProbDistro (const Target & t, double r) const;


private:
  class PairComposition: public std::pair<int,int>{
    bool Contains(int pdg) const noexcept{return pdg == first || pdg == second;}
    bool ContainsProton() const noexcept{return Contains(kPdgProton);}
    bool ContainsNeutron() const noexcept{return Contains(kPdgNeutron);}

  };

  std::vector<double> fContacts;
  std::vector<std::shared_ptr<Spline>> fUnivFunctions;
  std::vector<std::unique_ptr<Spline>> fInvCumulativeFunctions;
  std::vector<PairComposition> fPairComp;
  
  // options related to SRC pairs
  double fPCutOff;
  const NuclearModelI * fBaseModel;	
};

}         // genie namespace
#endif    // _LOCAL_FGM_H_

