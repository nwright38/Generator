//____________________________________________________________________________
/*!

\program testReadGNtuple

\brief   A simple test program to illustrate how to use the GENIE ROOT output
         TTrees of both formats.

         To run: testReadGNtuple -f filename
         where the filename points to a ROOT file containing a GENIE output
         TTree

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 02, 2005

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"

using std::string;
using namespace genie;

bool   testCommandLineArgs(int argc);
string getRootFilename    (int argc, char ** argv);
bool   checkRootFilename  (string filename);
void   readGHEPTrees      (TTree * tree); 

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments and get the ROOT filename
  string filename = getRootFilename(argc,argv);

  if ( !testCommandLineArgs(argc)   ) return 1;
  if ( !checkRootFilename(filename) ) return 2;

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(filename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("Main", pINFO) << "Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;

  //-- read & print the TTree entries [depends on the format]
  LOG("Main", pINFO) << "This ntuple's format is : "
                                   << NtpMCFormat::AsString(format);
  switch(format) {
     case kNFGHEP:
       // ntuple for passing data to other applications (eg the MC
       // simulation framework of an experiment)
       readGHEPTrees(tree);
       break;
     default:
       LOG("Main", pERROR) << "Cannot read this GENIE TTree format";
       return 3;
       break;
  }

  LOG("Main", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
bool testCommandLineArgs(int argc)
{
  if(argc!=3) {
   LOG("Main", pERROR) << "Not enough command line arguments";
   LOG("Main", pINFO)  << "Syntax: testReadGNtuple -f root_filename";
   return false;
  }
  return true;
}
//___________________________________________________________________
string getRootFilename(int argc, char ** argv)
{
  //-- Scan for filename from the command line argument (following -f)
  string filename="";
  for(int iarg = 0; iarg < argc-1; iarg++) {
   string argument(argv[iarg]);
   if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }
  return filename;
}
//___________________________________________________________________
bool checkRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("Main", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//___________________________________________________________________
void readGHEPTrees(TTree * tree)
{
// This ntuple contains a single TBranch with NtpMCEventRecord
// objects in its leaves

  if(!tree) return;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over TTree NtpMC records, get & print the events

  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("Main", pINFO) << rec_header;
    LOG("Main", pINFO) << event;

    mcrec->Clear();
  }
}
//___________________________________________________________________