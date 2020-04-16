//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file B4bEventAction.hh
/// \brief Definition of the B4bEventAction class

#ifndef B4bEventAction_h
#define B4bEventAction_h 1

#include "G4UserEventAction.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4SteppingManager.hh"

#include "globals.hh"

#include "G4Step.hh"

// for root tree
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TROOT.h"

#include <stdlib.h>     /* getenv */

#include "CaloDataStruc.h"

using namespace std;

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in B4bRunData object.

class B4bEventAction : public G4UserEventAction
{
  public:
    B4bEventAction(B4DetectorConstruction* det,B4PrimaryGeneratorAction* prim);
    virtual ~B4bEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AccumulateCaloHits(CaloStepData aHit);
    void FillSecondaries(const G4Step* step);
    void StepAnalysis(const G4Step* step);

  private:
    // methods
    B4DetectorConstruction* fDetector;
    B4PrimaryGeneratorAction* primary;

    void   clearTTreeVectors();   
    void   getCellSize();
    double TimeTOFAdjusted(double, double);

    // control parameters...
    string _CaloXG4OutName; //"caloxTree01.root"  # output root tree file name.
    int _CaloXG4RunNumber;  //=1   # run number
    int _CaloXG4EventNumber; //=1 # initial event number
    int _CaloXG4SeedA;       //=1       # 0=default (fixed), 1 auto, >1 manual seed setting
    int _CaloXG4SeedB;       //=2019    # second seed for manual seed setting (CaloXG4SeedA>1)
    int _CaloXG4SaveSeeds;   //=1   # 0=not save,  1=save seeds. 
    int _CaloXG4map1CellEdepON;    //=1     # Hit Map 0=off, 1=on
    int _CaloXG4map2PidEdepON;  //=1   # Edep per particle (Pid) 0=off, 1=on
    int _CaloXG4mapTcutMin;   //=100000 # time cut on Map creation  (psec)  
    double  d_CaloXG4mapTcutMin;   // above in nsec  
    int _CaloXG4mapTcutMax;   //=100000 # time cut on Map creation  (psec)  
    double  d_CaloXG4mapTcutMax;   // above in nsec  
    int _CaloXG4map1CelEcutMinSave; //=10  # energy cut on Edep per cell during tree output (kev).
    double d_CaloXG4map1CelEcutMinSave; //=10  # above in kev


    int caloNx, caloNy, caloNz;
    double caloDx,caloDy,caloDz;
    double caloXhalf, caloYhalf, caloZhalf;
    double worldZhalf;

    // variables for histograms...
    std::vector<double> edepSum;

    //  Hit Map (Edep in Cell)...
    std::map<int, double> caloHits;
    std::map<int, double> caloHits1;
    std::map<int, double> caloHits4;
    std::map<int, double> caloHits5;
    std::map<int, double> caloHits5em;
    std::map<int, double> caloHits10;
    std::map<int, double>::iterator caloHitsiter;

    //  Edep per Pid...
    std::map<int, double> tkPidEdep;
    std::map<int, double>::iterator tkPidEdepiter;
   
    //  Leakage kinetic energy per Pid
    std::map<int, double> leakPidKinE;
    std::map<int, double>::iterator leakPidKinEiter;

    std::map<int, double> backLeakPidKinE;
    std::map<int, double>::iterator backLeakPidKinEiter;

    std::map<int, double> frontLeakPidKinE;
    std::map<int, double>::iterator frontLeakPidKinEiter;

    // histograming and ntuple (ttree)...
      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

    TFile *fout;
    TTree *tree;

    int mRun;
    int mEvent;

    int mNxCell;     // number of cells in x
    int mNyCell;     // number of cells in y
    int mNzCell;     // number of cells in z
    float mDxCell;   // sell size in x (full length)
    float mDyCell;   // sell size in y (full length)
    float mDzCell;   // sell size in z (full length)
    float mzhalfWorld;  // half length of World Volume in Z

    std::vector<int>    mGenPID;   // pdg code
    std::vector<float>  mGenPx;    // MeV
    std::vector<float>  mGenPy;    // MeV
    std::vector<float>  mGenPz;    // MeV
    std::vector<float>  mGenE;    // MeV
    std::vector<float>  mGenMass;     // MeV

    std::vector<std::string> mVtxProcName ; // process name
    std::vector<int>    mVtxType;  // tbd
    std::vector<float>  mVtxX;     // mm
    std::vector<float>  mVtxY;     // mm
    std::vector<float>  mVtxZ;     // mm
    std::vector<float>  mVtxTime;  // ns
    std::vector<int>    mVtx2Part;  //  index of the first particles from this vertex
    std::vector<int>    mVtxNPart;  //  number of particles

    std::vector<int>    mPart2Vtx;  //  index of vertex which this particle associates with
    std::vector<int>    mPartPID;   // pdg code
    std::vector<float>  mPartPx;    // MeV
    std::vector<float>  mPartPy;    // MeV
    std::vector<float>  mPartPz;    // MeV
    std::vector<float>  mPartMass;  // MeV
    std::vector<float>  mPartKinE;    // Kinetic Energy (MeV)
    std::vector<float>  mPartTime;  // ns  (production time)

    std::vector<int>   mT0HitID;
    std::vector<float> mT0HitEdep;
    std::vector<float> mT0HitEsum;
    std::vector<int>   mT0HitCounts;

    std::vector<int>   mT1HitID;
    std::vector<float> mT1HitEdep;
    std::vector<float> mT1HitEsum;
    std::vector<int>   mT1HitCounts;

    std::vector<int>   mT4HitID;
    std::vector<float> mT4HitEdep;
    std::vector<float> mT4HitEsum;
    std::vector<int>   mT4HitCounts;

    std::vector<int>   mT10HitID;
    std::vector<float> mT10HitEdep;
    std::vector<float> mT10HitEsum;
    std::vector<int>   mT10HitCounts;

    std::vector<int>   mT5HitID;
    std::vector<float> mT5HitEdep;
    std::vector<float> mT5HitEsum;
    std::vector<int>   mT5HitCounts;

    std::vector<int>   mT5emHitID;
    std::vector<float> mT5emHitEdep;
    std::vector<float> mT5emHitEsum;
    std::vector<int>   mT5emHitCounts;

    std::vector<int>   mTk0Pid;
    std::vector<float> mTk0Edep;
    std::vector<float> mTk0Esum;
    std::vector<int>   mTk0Counts;

    std::vector<int>   mLeakPid;
    std::vector<float> mLeakKinE;
    std::vector<float> mLeakSum;
    std::vector<int>   mLeakCounts;

    std::vector<int>   mBackLeakPid;
    std::vector<float> mBackLeakKinE;
    std::vector<float> mBackLeakSum;
    std::vector<int>   mBackLeakCounts;

    std::vector<int>   mFrontLeakPid;
    std::vector<float> mFrontLeakKinE;
    std::vector<float> mFrontLeakSum;
    std::vector<int>   mFrontLeakCounts;

};


                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
inline void B4bEventAction::AccumulateCaloHits(CaloStepData aHit){

   std::cout<<"B4bEventAction::AccumulateCaloHits  track "<<aHit.trackid<<"  pid "<<aHit.pid<<"  edep "<<aHit.edep<<std::endl;
}
*/

#endif

    
