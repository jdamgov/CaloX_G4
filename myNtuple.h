//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 20 11:18:19 2019 by ROOT version 6.16/00
// from TTree tree/CaloX Tree
// found on file: caloxTree01_30GeV_10events.root
//////////////////////////////////////////////////////////

#ifndef myNtuple_h
#define myNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class myNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nGenPar;
   Float_t         nGenPt;
   vector<int>     *hitIndex;
   vector<float>   *hitEdep;

   // List of branches
   TBranch        *b_nGenPar;   //!
   TBranch        *b_nGenPt;   //!
   TBranch        *b_hitIndex;   //!
   TBranch        *b_hitEdep;   //!

   myNtuple(TTree *tree=0);
   virtual ~myNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myNtuple_cxx
myNtuple::myNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("caloxTree01_30GeV_10events.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("caloxTree01_30GeV_10events.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

myNtuple::~myNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hitIndex = 0;
   hitEdep = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
   fChain->SetBranchAddress("nGenPt", &nGenPt, &b_nGenPt);
   fChain->SetBranchAddress("hitIndex", &hitIndex, &b_hitIndex);
   fChain->SetBranchAddress("hitEdep", &hitEdep, &b_hitEdep);
   Notify();
}

Bool_t myNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myNtuple_cxx
