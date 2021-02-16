//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Author: Joe Johnston, University of Pittsburgh (Advisor Steven Dytman)

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/GeneralContactFormalism.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
GeneralContactFormalism::GeneralContactFormalism() :
NuclearModelI("genie::GeneralContactFormalism")
{

}
//____________________________________________________________________________
GeneralContactFormalism::GeneralContactFormalism(string config) :
NuclearModelI("genie::GeneralContactFormalism", config)
{

}
//____________________________________________________________________________
GeneralContactFormalism::~GeneralContactFormalism()
{

}
//____________________________________________________________________________
bool GeneralContactFormalism::GenerateNucleon(const Target & target,
				      double hitNucleonRadius) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  //-- set fermi momentum vector
  //
  TH1D * prob = this->ProbDistro(target,hitNucleonRadius);
  if(!prob) {
    LOG("GeneralContactFormalism", pNOTICE)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }
  double p = prob->GetRandom();
  delete prob;
  LOG("GeneralContactFormalism", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  fCurrMomentum.SetXYZ(px,py,pz);

  //-- set removal energy
  //
  int Z = target.Z();
  map<int,double>::const_iterator it = fNucRmvE.find(Z);
  if(it != fNucRmvE.end()) fCurrRemovalEnergy = it->second;
  else fCurrRemovalEnergy = nuclear::BindEnergyPerNucleon(target);

  return true;
}
//____________________________________________________________________________
double GeneralContactFormalism::Prob(double p, double w, const Target & target,
			     double hitNucleonRadius) const
{
  if(w<0) {
    TH1D * prob = this->ProbDistro(target, hitNucleonRadius);
    int bin = prob->FindBin(p);
    double y  = prob->GetBinContent(bin);
    double dx = prob->GetBinWidth(bin);
    double pr  = y*dx;
    delete prob;
    return pr;
  }
  return 1;
}
//____________________________________________________________________________
// *** The TH1D object must be deleted after it is used ***
TH1D * GeneralContactFormalism::ProbDistro(const Target & target, double r) const
{
  LOG("GeneralContactFormalism", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString()
	     << ", Nucleon Radius = " << r;
  LOG("GeneralContactFormalism", pNOTICE)
             << ", P(max) = " << fPMax;

  assert(target.HitNucIsSet());

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));

  // bool is_p = pdg::IsProton(nucleon_pdgc);
  // double numNuc = (is_p) ? (double)target.Z():(double)target.N();

  // Calculate Fermi Momentum using Local FG equations
  double KF = LocalFermiMomentum( target, nucleon_pdgc, r ) ; 

  LOG("GeneralContactFormalism",pNOTICE) << "KF = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;

  // Do not include nucleon correlation tail
  //double R  = 1. / (1.- KF/4.);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GeneralContactFormalism", pDEBUG) << "a  = " << a;
  LOG("GeneralContactFormalism", pDEBUG) << "C  = " << C;
  //LOG("GeneralContactFormalism", pDEBUG) << "R  = " << R;
#endif

  //-- create the probability distribution

  int npbins = (int) (1000*fPMax);
  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);






  double dp = fPMax / (npbins-1);
//  double iC = (C>0) ? 1./C : 0.; // unused variables
//  double kfa_pi_2 = TMath::Power(KF*a/kPi,2); // unused variables

  for(int i = 0; i < npbins; i++) {
     double p  = i * dp;
     double p2 = TMath::Power(p,2);

     // use expression with fSRC_Fraction to allow the possibility of 
     // using the Correlated Fermi Gas Model with a high momentum tail

     // calculate |phi(p)|^2
     double phi2 = 0;
        if (p <= KF){
            phi2 = (1./(4*kPi)) * (3/TMath::Power(KF,3.)) * ( 1 - fSRC_Fraction );
        }else if( p > KF && p < fPCutOff ){
            phi2 = (1./(4*kPi)) * ( fSRC_Fraction / (1./KF - 1./fPCutOff) ) / TMath::Power(p,4.);
        }

     // calculate probability density : dProbability/dp
     double dP_dp = 4*kPi * p2 * phi2;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GeneralContactFormalism", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob;
}
//____________________________________________________________________________
double GeneralContactFormalism::LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const {

  assert(pdg::IsProton(nucleon_pdg) || pdg::IsNeutron(nucleon_pdg)) ;

  bool is_p = pdg::IsProton(nucleon_pdg);
  double numNuc = (double) ( (is_p) ? t.Z() : t.N() );

  //  double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
  
  double kF = TMath::Power( 3*kPi2*numNuc*genie::utils::nuclear::Density( radius, t.A() ),
			   1.0/3.0 ) 
    / genie::units::fermi ;

  return kF ;
}
//____________________________________________________________________________
void GeneralContactFormalism::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GeneralContactFormalism::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void GeneralContactFormalism::LoadConfig(void)
{

  NuclearModelI::LoadConfig() ;

  double contact = 0;
  this->GetParamDef("Contact-pn1", contact, 0.);

  string filename;
  this->GetParam("UniversalFunction-pn1", filename);
  if(filename.find('/') == string::npos){ 
    filename = string(gSystem->Getenv("GENIE")) + "/data/evgen/nucl/contacts/" + filename;
  }

  Spline s(filename);

  std::vector<double> x,y;
  for(unsigned_int i = 0; i < s.NKnots(); ++i){
    x.push_back(s.GetKnotX(i)/units::fm);
    y.push_back(s.GetKnotY(i)*units::fm3*4*constants::kPi);

  }

  fUnivFunctions.push_back(new Spline(x.size(),x.data(),y.data()));

  
//Evaluate inv cumulative fn 
//Convert to GeV and to radial distribution
//Normalize to 1
//  do the integral

  for(unsigned int i = 0; i < fUnivFunctions.size(); ++i){
    
    GSLIntegrator integrator(kADAPTIVE);
    GCF::PDF_Pr(fUnivFunctions[i]);
    double integral = integrator.Integral(PDF_Pr,xmin,fUnivFunctions[i]->Xmax());
    fUnivFunctions[i]->Divide(integral);  
    //define inv cumulative fn
    //define arrays x and y,  use for a spline, (make last point 1?) 

 }
  
 


  this->GetParamDef("LFG-MomentumMax", fPMax, 1.0);
  assert(fPMax > 0);

  this->GetParamDef("SRC-Fraction", fSRC_Fraction, 0.0);
  this->GetParam("LFG-MomentumCutOff", fPCutOff);

  if (fPCutOff > fPMax) {
        LOG("GeneralContactFormalism", pFATAL) << "Momentum CutOff greater than Momentum Max";
        exit(78);
  }


  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
  //

  for(int Z=1; Z<140; Z++) {
    for(int A=Z; A<3*Z; A++) {
      ostringstream key ;
      int pdgc = pdg::IonPdgCode(A,Z);
      key << "RFG-NucRemovalE@Pdg=" << pdgc;
      RgKey rgkey = key.str();
      double eb ;
      if ( GetParam( rgkey, eb, false ) ) {
    	eb = TMath::Max(eb, 0.);
        LOG("GeneralContactFormalism", pINFO) << "Nucleus: " << pdgc << " -> using Eb =  " << eb << " GeV";
        fNucRmvE.insert(map<int,double>::value_type(Z,eb));
      }
    }
  }
}
//____________________________________________________________________________
