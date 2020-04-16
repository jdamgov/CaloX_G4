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
/// \file B4bEventAction.cc
/// \brief Implementation of the B4bEventAction class

#include "B4bEventAction.hh"
#include "B4bRunData.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "G4Step.hh"

#include "B4PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include <iomanip>

// -- for root ntuple --
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TText.h"
#include "TPaveText.h"
#include "TROOT.h"

// -- for CaloX data --
#include "CaloDataStruc.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::B4bEventAction(B4DetectorConstruction* det,B4PrimaryGeneratorAction* prim)
 : G4UserEventAction(),fDetector(det),primary(prim)
{  

   cout<<"initializing CaloXTree..."<<endl;

   //  get control parameters from linux env varaiabes...
   _CaloXG4OutName="caloxTree01.root";
     char* _param1;
     _param1 = getenv("CaloXG4OutName");
     if(_param1 != NULL) {string ss(_param1); _CaloXG4OutName=ss;}
   
   _CaloXG4RunNumber=1;
     char* _param2;
     _param2 = getenv("CaloXG4RunNumber");
     if(_param2 != NULL) {_CaloXG4RunNumber=atoi(_param2);}

   _CaloXG4EventNumber=1;
     char* _param3;
     _param3 = getenv("CaloXG4EventNumber");
     if(_param3 != NULL) {_CaloXG4EventNumber=atoi(_param3);}

////////
   _CaloXG4SeedA=1;
     char* _param4;
     _param4 = getenv("CaloXG4SeedA");
     if(_param4 != NULL) {_CaloXG4SeedA=atoi(_param4);}

   _CaloXG4SeedB=1;
     char* _param5;
     _param5 = getenv("CaloXG4SeedB");
     if(_param5 != NULL) {_CaloXG4SeedB=atoi(_param5);}

   _CaloXG4SaveSeeds=1;
     char* _param6;
     _param6 = getenv("CaloXG4SaveSeeds");
     if(_param6 != NULL) {_CaloXG4SaveSeeds=atoi(_param6);}

   _CaloXG4map1CellEdepON=1;
     char* _param7;
     _param7 = getenv("CaloXG4map1CellEdepON");
     if(_param7 != NULL) {_CaloXG4map1CellEdepON=atoi(_param7);}

   _CaloXG4map2PidEdepON=1;
     char* _param8;
     _param8 = getenv("CaloXG4map2PidEdepON");
     if(_param8 != NULL) {_CaloXG4map2PidEdepON=atoi(_param8);}

   _CaloXG4mapTcutMin=0 ;  // in psec
     char* _param9;
     _param9 = getenv("CaloXG4mapTcutMin");
     if(_param9 != NULL) {_CaloXG4mapTcutMin=atoi(_param9);
                d_CaloXG4mapTcutMin=double(_CaloXG4mapTcutMin)/1000.0;}  // in nsec 

   _CaloXG4mapTcutMax=100000;  // in psec  (default 100 ns)
     char* _param9a;
     _param9a = getenv("CaloXG4mapTcutMax");
     if(_param9a != NULL) {_CaloXG4mapTcutMax=atoi(_param9a);
                d_CaloXG4mapTcutMax=double(_CaloXG4mapTcutMax)/1000.0;}  // in nsec (internal) 

   _CaloXG4map1CelEcutMinSave=10;  // in kev
     char* _param10;
     _param10 = getenv("CaloXG4map1CelEcutMinSave");
     if(_param10 != NULL) {_CaloXG4map1CelEcutMinSave=atoi(_param10);
                d_CaloXG4map1CelEcutMinSave=_CaloXG4map1CelEcutMinSave/1000.0;} // in MeV (intenal)
   //
   //   random number seeds...
   //
   if(_CaloXG4SeedA>0) {
      if(_CaloXG4SeedA==1) {
         long seeds[2];
         time_t systime = time(NULL);
         seeds[0] = (long) (systime+_CaloXG4RunNumber);
         seeds[1] = (long) (systime*G4UniformRand());
         G4Random::setTheSeeds(seeds);
         G4Random::showEngineStatus();
      }
      if(_CaloXG4SeedA>1) {
         long seeds[2];
         seeds[0] = (long) _CaloXG4SeedA;
         seeds[1] = (long) _CaloXG4SeedB;
         G4Random::setTheSeeds(seeds);
         G4Random::showEngineStatus();
      }
   }

   if(_CaloXG4SaveSeeds>0) {
       string rndmFile="BeginOfRun_"+to_string(_CaloXG4RunNumber)+".rndm";
       G4Random::saveEngineStatus(rndmFile.c_str());
   }
   //
   //   define output root tree file...
   //
   // string outname="caloxTree01.root";
   fout=new TFile(_CaloXG4OutName.c_str(),"recreate");

   histo1D["edepALL"]=new TH1D("edepAll","Total Edep (MeV)",1000,0.,500000.);
   histo1D["edepT5"]=new TH1D("edepT5","Total Edep (MeV), T<5ns",1000,0.,500000.);
   histo1D["edepT10"]=new TH1D("edepT10","Total Edep (MeV), T<10ns",1000,0.,500000.);
   histo1D["edepT50"]=new TH1D("edepT50","Total Edep (MeV), T<50ns",1000,0.,500000.);
   histo1D["edepT100"]=new TH1D("edepT100","Total Edep (MeV), T<100ns",1000,0.,500000.);
   histo1D["edepTime"]=new TH1D("edepTime","Time of Edep (ns)",1000,0.,100.);
   histo1D["edepTimeEnWt"]=new TH1D("edepTimeEnWt","Time of Edep (ns) (edep weighted)",1000,0.,100.);

   std::vector<double> v(10,0.0);
   edepSum=v;

      // define root tree...

   tree=new TTree("tree","CaloX Tree");

   tree->Branch("Run" , &mRun);
   tree->Branch("Event", &mEvent);

   tree->Branch("NxCell" , &mNxCell);
   tree->Branch("NyCell" , &mNyCell);
   tree->Branch("NzCell" , &mNzCell);
   tree->Branch("DxCell" , &mDxCell);
   tree->Branch("DyCell" , &mDyCell);
   tree->Branch("DzCell" , &mDzCell);
   tree->Branch("zhalfWorld" , &mzhalfWorld);


   tree->Branch("GenPID"   , &mGenPID);
   tree->Branch("GenPx"    , &mGenPx);
   tree->Branch("GenPy"    , &mGenPy);
   tree->Branch("GenPz"    , &mGenPz);
   tree->Branch("GenE"     , &mGenE);
   tree->Branch("GenMass"  , &mGenMass);

   tree->Branch("VtxProcName", &mVtxProcName);
   tree->Branch("VtxType"    , &mVtxType);
   tree->Branch("VtxX"       , &mVtxX);
   tree->Branch("VtxY"       , &mVtxY);
   tree->Branch("VtxZ"       , &mVtxZ);
   tree->Branch("VtxTime"    , &mVtxTime);
   tree->Branch("VtxPartPtr" , &mVtx2Part);
   tree->Branch("VtxNPart"   , &mVtxNPart);

   tree->Branch("PartVtxPtr"  , &mPart2Vtx);
   tree->Branch("PartPID"  , &mPartPID);
   tree->Branch("PartPx"  , &mPartPx);
   tree->Branch("PartPy"  , &mPartPy);
   tree->Branch("PartPz"  , &mPartPz);
   tree->Branch("PartMass"  , &mPartMass);
   tree->Branch("PartKinE"  , &mPartKinE);
   tree->Branch("PartTime"  , &mPartTime);


   tree->Branch("T0HitID" , &mT0HitID);
   tree->Branch("T0HitEdep"  , &mT0HitEdep);
   tree->Branch("T0HitEsum"  , &mT0HitEsum);
   tree->Branch("T0HitCounts"  , &mT0HitCounts);

   tree->Branch("T1HitID" , &mT1HitID);
   tree->Branch("T1HitEdep"  , &mT1HitEdep);
   tree->Branch("T1HitEsum"  , &mT1HitEsum);
   tree->Branch("T1HitCounts"  , &mT1HitCounts);

   tree->Branch("T4HitID" , &mT4HitID);
   tree->Branch("T4HitEdep"  , &mT4HitEdep);
   tree->Branch("T4HitEsum"  , &mT4HitEsum);
   tree->Branch("T4HitCounts"  , &mT4HitCounts);

   tree->Branch("T5HitID" , &mT5HitID);
   tree->Branch("T5HitEdep"  , &mT5HitEdep);
   tree->Branch("T5HitEsum"  , &mT5HitEsum);
   tree->Branch("T5HitCounts"  , &mT5HitCounts);

   tree->Branch("T5emHitID" , &mT5emHitID);
   tree->Branch("T5emHitEdep"  , &mT5emHitEdep);
   tree->Branch("T5emHitEsum"  , &mT5emHitEsum);
   tree->Branch("T5emHitCounts"  , &mT5emHitCounts);

   tree->Branch("T10HitID" , &mT10HitID);
   tree->Branch("T10HitEdep"  , &mT10HitEdep);
   tree->Branch("T10HitEsum"  , &mT10HitEsum);
   tree->Branch("T10HitCounts"  , &mT10HitCounts);

   tree->Branch("Tk0Pid" , &mTk0Pid);
   tree->Branch("Tk0Edep" , &mTk0Edep);
   tree->Branch("Tk0Esum" , &mTk0Esum);
   tree->Branch("Tk0Counts" , &mTk0Counts);

   tree->Branch("LeakPid" , &mLeakPid);
   tree->Branch("LeakKinE" , &mLeakKinE);
   tree->Branch("LeakSum" , &mLeakSum);
   tree->Branch("LeakCounts" , &mLeakCounts);

   tree->Branch("BackLeakPid" , &mBackLeakPid);
   tree->Branch("BackLeakKinE" , &mBackLeakKinE);
   tree->Branch("BackLeakSum" , &mBackLeakSum);
   tree->Branch("BackLeakCounts" , &mBackLeakCounts);

   tree->Branch("FrontLeakPid" , &mFrontLeakPid);
   tree->Branch("FrontLeakKinE" , &mFrontLeakKinE);
   tree->Branch("FrontLeakSum"  , &mFrontLeakSum);
   tree->Branch("FrontLeakCounts" , &mFrontLeakCounts);
   // set event counter.
   mEvent=0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::~B4bEventAction()
{
  fout->Write();
  fout->Close();

  if(_CaloXG4SaveSeeds>0) {
     G4Random::showEngineStatus();
     string rndmFile="endOfRun_"+to_string(_CaloXG4RunNumber)+".rndm";
     G4Random::saveEngineStatus(rndmFile.c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  clearTTreeVectors();
  getCellSize();

  if(caloHits.size()>0) caloHits.clear();
  if(caloHits1.size()>0) caloHits1.clear();
  if(caloHits4.size()>0) caloHits4.clear();
  if(caloHits5.size()>0) caloHits5.clear();
  if(caloHits5em.size()>0) caloHits5em.clear();
  if(caloHits10.size()>0) caloHits10.clear();
  if(tkPidEdep.size()>0) tkPidEdep.clear();
  if(leakPidKinE.size()>0) leakPidKinE.clear();
  if(backLeakPidKinE.size()>0) backLeakPidKinE.clear();
  if(frontLeakPidKinE.size()>0) frontLeakPidKinE.clear();

  // clear some vectors...
  std::fill(edepSum.begin(),edepSum.end(),0.0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::EndOfEventAction(const G4Event* event)
{   
   // std::cout<<"EndOfEventAction... edepSum[0]  "<<edepSum[0]<<std::endl;
   histo1D["edepALL"]->Fill(edepSum[0]);
   histo1D["edepT5"]->Fill(edepSum[1]);
   histo1D["edepT10"]->Fill(edepSum[2]);
   histo1D["edepT50"]->Fill(edepSum[3]);
   histo1D["edepT100"]->Fill(edepSum[4]);

   mRun=_CaloXG4RunNumber;
   mEvent++;

   // std::cout<<"B4bEventAction::EndOfEventAction...   loop over hits "<<std::endl;

   // double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
   mGenPID.push_back(primary->GetParticleGun()->GetParticleDefinition()->GetPDGEncoding());
   mGenPx.push_back(primary->GetParticleGun()->GetParticleMomentumDirection().x());
   mGenPy.push_back(primary->GetParticleGun()->GetParticleMomentumDirection().y());
   mGenPz.push_back(primary->GetParticleGun()->GetParticleMomentumDirection().z());
   mGenE.push_back(primary->GetParticleGun()->GetParticleEnergy());
   mGenMass.push_back(primary->GetParticleGun()->GetParticleDefinition()->GetPDGMass());

   // int i=0;
   double esum=0.0;
   for (caloHitsiter=caloHits.begin(); caloHitsiter != caloHits.end(); caloHitsiter++) {
     // std::cout<<i<<"   "<<caloHitsiter->first<<"        "<<caloHitsiter->second<<std::endl;
     esum=esum+caloHitsiter->second;
    //  i++;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT0HitID.push_back(caloHitsiter->first);
        mT0HitEdep.push_back(caloHitsiter->second);
     }
   }
   mT0HitEsum.push_back(esum);
   mT0HitCounts.push_back(mT0HitID.size());

   esum=0.0;
   for (caloHitsiter=caloHits5.begin(); caloHitsiter != caloHits5.end(); caloHitsiter++) {
     esum=esum+caloHitsiter->second;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT5HitID.push_back(caloHitsiter->first);
        mT5HitEdep.push_back(caloHitsiter->second);
     }
   }
   mT5HitEsum.push_back(esum);
   mT5HitCounts.push_back(mT5HitID.size());

   esum=0.0;
   for (caloHitsiter=caloHits5em.begin(); caloHitsiter != caloHits5em.end(); caloHitsiter++) {
     esum=esum+caloHitsiter->second;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT5emHitID.push_back(caloHitsiter->first);
        mT5emHitEdep.push_back(caloHitsiter->second);
     }
   }
   mT5emHitEsum.push_back(esum);
   mT5emHitCounts.push_back(mT5emHitID.size());

   esum=0.0;
   for (caloHitsiter=caloHits10.begin(); caloHitsiter != caloHits10.end(); caloHitsiter++) {
     esum=esum+caloHitsiter->second;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT10HitID.push_back(caloHitsiter->first);
        mT10HitEdep.push_back(caloHitsiter->second);
     }
   }
   mT10HitEsum.push_back(esum);
   mT10HitCounts.push_back(mT10HitID.size());

   esum=0.0;
   for (caloHitsiter=caloHits1.begin(); caloHitsiter != caloHits1.end(); caloHitsiter++) {
     esum=esum+caloHitsiter->second;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT1HitID.push_back(caloHitsiter->first);
        mT1HitEdep.push_back(caloHitsiter->second);
     }
   }
   mT1HitEsum.push_back(esum);
   mT1HitCounts.push_back(mT1HitID.size());

   esum=0.0;
   for (caloHitsiter=caloHits4.begin(); caloHitsiter != caloHits4.end(); caloHitsiter++) {
     esum=esum+caloHitsiter->second;
     if(caloHitsiter->second>d_CaloXG4map1CelEcutMinSave) {
        mT4HitID.push_back(caloHitsiter->first);
        mT4HitEdep.push_back(caloHitsiter->second);
     }
   }
   mT4HitEsum.push_back(esum);
   mT4HitCounts.push_back(mT4HitID.size());


   // std::cout<<"total energy deposit = "<<esum<<"  in "<<caloHits.size()<<" hits."<<std::endl;
   // std::cout<<"   mHitID    size  "<<mT0HitID.size()<<std::endl;
   // std::cout<<"   mHitEdep  size  "<<mT0HitEdep.size()<<std::endl;

   //   map2:  Pid vs Edep
   esum=0.0;
   for (tkPidEdepiter=tkPidEdep.begin(); tkPidEdepiter != tkPidEdep.end(); tkPidEdepiter++) {
     mTk0Pid.push_back(tkPidEdepiter->first);
     mTk0Edep.push_back(tkPidEdepiter->second);
     esum=esum+tkPidEdepiter->second;
  }
  mT0HitEsum.push_back(esum);
  mTk0Counts.push_back(mTk0Pid.size());

  //    energy leakage, all, back, front
  esum=0.0;
  for (leakPidKinEiter=leakPidKinE.begin(); leakPidKinEiter != leakPidKinE.end(); leakPidKinEiter++) {
     mLeakPid.push_back(leakPidKinEiter->first);
     mLeakKinE.push_back(leakPidKinEiter->second);
     esum=esum+leakPidKinEiter->second;
  }
  mLeakSum.push_back(esum);
  mLeakCounts.push_back(mLeakPid.size());

  esum=0.0;
  for (backLeakPidKinEiter=backLeakPidKinE.begin(); backLeakPidKinEiter != backLeakPidKinE.end(); backLeakPidKinEiter++) {
     mBackLeakPid.push_back(backLeakPidKinEiter->first);
     mBackLeakKinE.push_back(backLeakPidKinEiter->second);
     esum=esum+backLeakPidKinEiter->second;
  }
  mBackLeakSum.push_back(esum);
  mBackLeakCounts.push_back(mBackLeakPid.size());

  esum=0.0;
  for (frontLeakPidKinEiter=frontLeakPidKinE.begin(); frontLeakPidKinEiter != frontLeakPidKinE.end(); frontLeakPidKinEiter++) {
     mFrontLeakPid.push_back(frontLeakPidKinEiter->first);
     mFrontLeakKinE.push_back(frontLeakPidKinEiter->second);
     esum=esum+frontLeakPidKinEiter->second;
  }
  mFrontLeakSum.push_back(esum);
  mFrontLeakCounts.push_back(mBackLeakPid.size());

   tree->Fill();

} //  end of B4bEventAction::EndOfEventAction  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::AccumulateCaloHits(CaloStepData aHit){
      // std::cout<<"B4bEventAction::AccumulateCaloHits  track "<<aHit.trackid<<"  pid "<<aHit.pid<<"  edep "<<aHit.edep<<std::endl;
   // if(aHit.edep>0.000001) { 
   if(aHit.edep>0.0) { 
      double tA=TimeTOFAdjusted(aHit.globaltime,aHit.z);
      if(tA>d_CaloXG4mapTcutMin && tA<d_CaloXG4mapTcutMax) {
         int kx=(aHit.x+caloXhalf)/caloDx;
         int ky=(aHit.y+caloYhalf)/caloDy;
         int kz=(aHit.z+caloZhalf)/caloDz;
         int index=kx*1000000+ky*1000+kz;
         if(_CaloXG4map1CellEdepON>0) { 
            caloHits[index]=caloHits[index]+aHit.edep;
            if( tA <=1.) caloHits1[index]=caloHits1[index]+aHit.edep;
            if( tA <=5.) caloHits5[index]=caloHits5[index]+aHit.edep;
            if( tA <=5. && abs(aHit.pid)==11) caloHits5em[index]=caloHits5em[index]+aHit.edep;
            if( tA>1 && tA <=4.0) caloHits4[index]=caloHits4[index]+aHit.edep;
            if( tA>4.0 && tA <=10.0) caloHits10[index]=caloHits10[index]+aHit.edep;
         }
         //
         if(_CaloXG4map2PidEdepON>0) {
            tkPidEdep[aHit.pid]=tkPidEdep[aHit.pid]+aHit.edep;
         }
      } // end of if(tA>d_CaloXG4mapTcutMin && tA<d_CaloXG4mapTcutMax)
   }  // end of if(aHit.edep>0.0)

}

// -----------------------------------------------------------------------
void B4bEventAction::FillSecondaries(const G4Step* step){
      // THis is called from SteppingAction if(tkstatus==fStopAndKill && absPdgCode>100)
      // futher selection of interactions should be done in this code.

      // Two options to get secondaries- (1) only secondaries produced at the 
      // current step,  (2) all secondaries includeing delta-ray electrons, etc.
      // particles produced along the track up to the current point.
      const std::vector<const G4Track*>* fSecondary=step->GetSecondaryInCurrentStep();
      // const G4TrackVector* fSecondary=step->GetSecondary(); // all seondaries

      G4Track* track = step->GetTrack();
      G4TrackStatus tkstatus=track->GetTrackStatus();
      G4String processName=step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

      const G4DynamicParticle* dynamicParticle= track ->GetDynamicParticle();
      G4int pdgcode=dynamicParticle->GetPDGcode();
      // G4int absPdgCode=abs(pdgcode);
      G4ParticleDefinition* particle = dynamicParticle->GetDefinition();
      G4String particleName= particle ->GetParticleName();
      G4ThreeVector posB = step->GetPostStepPoint()->GetPosition();

      //

      bool fillFlag=false;  // save one enegetic interactions...
      if(step->GetPreStepPoint()->GetKineticEnergy()>50.0 ) fillFlag=true; 
      if(fillFlag) {
         int vertexId=mVtxType.size();
         mVtxProcName.push_back(processName);
         mVtxType.push_back(0);
         mVtxX.push_back(posB.x());    
         mVtxY.push_back(posB.y());    
         mVtxZ.push_back(posB.z());    
         double t=TimeTOFAdjusted(track->GetGlobalTime(),posB.z());
         if(t>999.0) t=999.0;
         mVtxTime.push_back(t);    
         mVtx2Part.push_back(mPart2Vtx.size());
         int n=0;
         //  fill current track, at production point
                mPart2Vtx.push_back(vertexId);
                mPartPID.push_back(pdgcode);
                mPartPx.push_back(track->GetVertexMomentumDirection().x());
                mPartPy.push_back(track->GetVertexMomentumDirection().y());
                mPartPz.push_back(track->GetVertexMomentumDirection().z());
                mPartMass.push_back(track->GetDynamicParticle()->GetMass());
                mPartKinE.push_back(track->GetVertexKineticEnergy());
                // time when the parent particle was produced.
                double time1a=TimeTOFAdjusted(track->GetGlobalTime()-track->GetLocalTime()
                 ,track->GetVertexPosition().z());
                if(time1a>999.0) time1a=999.0;
                mPartTime.push_back(time1a);
         n++;
         //  fill current track, at previou step.
                mPart2Vtx.push_back(vertexId);
                mPartPID.push_back(pdgcode);
                mPartPx.push_back(step->GetPreStepPoint()->GetMomentumDirection().x());
                mPartPy.push_back(step->GetPreStepPoint()->GetMomentumDirection().y());
                mPartPz.push_back(step->GetPreStepPoint()->GetMomentumDirection().z());
                mPartMass.push_back(track->GetDynamicParticle()->GetMass());
                mPartKinE.push_back(step->GetPreStepPoint()->GetKineticEnergy());
                double time1b=TimeTOFAdjusted(track->GetGlobalTime(),track->GetPosition().z());
                if(time1b>999.0) time1b=999.0;
                mPartTime.push_back(time1b);
         n++;
         if(fSecondary->size()>0) {
            for (int lp1=0; lp1<fSecondary->size(); lp1++) {
                const G4Track* tk2= (*fSecondary)[lp1];
                G4int pdg2=tk2->GetDynamicParticle()->GetPDGcode();
                mPart2Vtx.push_back(vertexId);
                mPartPID.push_back(pdg2);
                mPartPx.push_back(tk2->GetMomentum().x());
                mPartPy.push_back(tk2->GetMomentum().y());
                mPartPz.push_back(tk2->GetMomentum().z());
                mPartMass.push_back(tk2->GetDynamicParticle()->GetMass());
                mPartKinE.push_back(tk2->GetKineticEnergy());
                double time2=TimeTOFAdjusted(tk2->GetGlobalTime(),tk2->GetPosition().z());
                if(time2>999.0) time2=999.0;
                mPartTime.push_back(time2);
                n++;
           }
         }
         mVtxNPart.push_back(n);   // store the number of secondaries in vtx block.
      }

      bool printFlag=false;
      // if(track->GetTrackID()==1) printFlag=true;
     if(printFlag) {
      std::cout<<"FillSecondaries.  processName="<<processName<<std::endl;
      std::cout<<" trackID "<<track->GetTrackID();
      std::cout<<"  "<<particleName;
      std::cout<<"   status "<<tkstatus;
      std::cout<<"  length "<<track->GetTrackLength();
      std::cout<<"  time "<<track->GetGlobalTime();
      std::cout<<"  proc="<<processName;
      std::cout<<"   Secondaries: "<<fSecondary->size();
      std::cout<<" "<<std::endl;

      if(fSecondary->size()>0) {
         std::cout << "    ++List of secondaries generated "
                << "(x,y,z,kE,t,PID):"
                << "  No. of secodaries = "
                << (*fSecondary).size() << std::endl;
         for (int lp1=0; lp1<fSecondary->size(); lp1++) {
                std::cout<<std::setprecision(3)
                << std::setw( 9)
                << (*fSecondary)[lp1]->GetPosition().x() << " x mm  "
                << std::setw( 9)
                << (*fSecondary)[lp1]->GetPosition().y() << " y mm  "
                << std::setw( 9)
                << (*fSecondary)[lp1]->GetPosition().z() << " z mm  "
                << std::setw( 9)
                << (*fSecondary)[lp1]->GetKineticEnergy() << " MeV "
                << std::setw( 9)
                << (*fSecondary)[lp1]->GetGlobalTime() << " ns "
                << std::setw(18)
                << (*fSecondary)[lp1]->GetDefinition()->GetParticleName() << std::endl;
        }
        // for checking output vector...
        for(int i=0; i<mPart2Vtx.size(); i++) {
          std::cout<<std::setprecision(3)
          <<std::setw(3)<<i
          <<std::setw(9)<<mPart2Vtx[i]<<"   "
          <<std::setw(9)<<mPartPID[i]<<" id "
          <<std::setw(9)<<mPartPx[i]<<" px  "
          <<std::setw(9)<<mPartPy[i]<<" py  "
          <<std::setw(9)<<mPartPz[i]<<" pz  "
          <<std::setw(9)<<mPartMass[i]<<" E   "
          <<std::setw(9)<<mPartKinE[i]<<" KinE "
          <<std::setw(9)<<mPartTime[i]<<" ns  "
          <<std::endl;
        }
   
     }
    }  // end of if(printFlag) 

}  // end of B4bEventAction::FillSecondaries.

// -----------------------------------------------------------------------
double B4bEventAction::TimeTOFAdjusted(double t, double  z){
    // return the time with subtraction of TOF along Z. 
    return t-(z+worldZhalf)/300.0;  // 300 mm/1 ns
}

// -----------------------------------------------------------------------
void B4bEventAction::StepAnalysis(const G4Step* step){
   double edep = step->GetTotalEnergyDeposit();
   G4Track* track = step->GetTrack();
   if(edep>0.0) {
      // G4Track* track = step->GetTrack();
      double tA=TimeTOFAdjusted(track->GetGlobalTime(),track->GetPosition().z());
      histo1D["edepTime"]->Fill(tA);
      histo1D["edepTimeEnWt"]->Fill(tA,edep);

      edepSum[0]=edepSum[0]+edep;
      if(tA<5.0) edepSum[1]=edepSum[1]+edep;
      if(tA<10.0) edepSum[2]=edepSum[2]+edep;
      if(tA<50.0) edepSum[3]=edepSum[3]+edep;
      if(tA<100.0) edepSum[4]=edepSum[4]+edep;
   }

  //  energy leakage...
  auto aPreStepPoint = step->GetPreStepPoint();
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto thisPhysical = touchable->GetVolume(); // mother
  auto thisName = thisPhysical->GetName();
  if(thisName.compare(0,5,"World")==0)  {
     G4bool just_entered = aPreStepPoint->GetStepStatus() == fGeomBoundary;
     if(just_entered) {
        const G4DynamicParticle* dynamicParticle= track ->GetDynamicParticle();
        G4int pdgcode=dynamicParticle->GetPDGcode();
        G4double kinEnergy=dynamicParticle->GetKineticEnergy();
        G4ThreeVector posA = step->GetPreStepPoint()->GetPosition();
        // cout<<"  just_entered "<<just_entered;
        // cout<<" pre(z) "<<posA.z();
        // cout<<" post(z) "<< track->GetPosition().z();
        // cout<<" caloZhalf "<<caloZhalf;
        // cout<<"  pdgcode "<< pdgcode<<endl;
        leakPidKinE[pdgcode]=leakPidKinE[pdgcode]+kinEnergy;
        if(posA.z()>(caloZhalf-0.1)) {
          backLeakPidKinE[pdgcode]=backLeakPidKinE[pdgcode]+kinEnergy;
        }
        if(posA.z()<(-caloZhalf+0.1)) {
          frontLeakPidKinE[pdgcode]=frontLeakPidKinE[pdgcode]+kinEnergy;
        } 
     }
  } 

}

// -----------------------------------------------------------------------
void B4bEventAction::clearTTreeVectors(){
   mNxCell=0;
   mNyCell=0;
   mNzCell=0;
   mDxCell=0.0;
   mDyCell=0.0;
   mDzCell=0.0;

   mGenPID.clear();
   mGenPx.clear();
   mGenPy.clear();
   mGenPz.clear();
   mGenE.clear();
   mGenMass.clear();

   mVtxProcName.clear();
   mVtxType.clear();
   mVtxX.clear();
   mVtxY.clear();
   mVtxZ.clear();
   mVtxTime.clear();
   mVtx2Part.clear();
   mVtxNPart.clear();

   mPart2Vtx.clear();
   mPartPID.clear();
   mPartPx.clear();
   mPartPy.clear();
   mPartPz.clear();
   mPartMass.clear();
   mPartKinE.clear();
   mPartTime.clear();

   mT0HitID.clear();
   mT0HitEdep.clear();
   mT0HitEsum.clear();
   mT0HitCounts.clear();

   mT10HitID.clear();
   mT10HitEdep.clear();
   mT10HitEsum.clear();
   mT10HitCounts.clear();

   mT5HitID.clear();
   mT5HitEdep.clear();
   mT5HitEsum.clear();
   mT5HitCounts.clear();

   mT5emHitID.clear();
   mT5emHitEdep.clear();
   mT5emHitEsum.clear();
   mT5emHitCounts.clear();

   mT1HitID.clear();
   mT1HitEdep.clear();
   mT1HitEsum.clear();
   mT1HitCounts.clear();

   mT4HitID.clear();
   mT4HitEdep.clear();
   mT4HitEsum.clear();
   mT4HitCounts.clear();

   mTk0Pid.clear();
   mTk0Edep.clear();
   mTk0Esum.clear();
   mTk0Counts.clear();

   mLeakPid.clear();
   mLeakKinE.clear();
   mLeakSum.clear();
   mLeakCounts.clear();

   mBackLeakPid.clear();
   mBackLeakKinE.clear();
   mBackLeakSum.clear();
   mBackLeakCounts.clear();

   mFrontLeakPid.clear();
   mFrontLeakKinE.clear();
   mFrontLeakSum.clear();
   mFrontLeakCounts.clear();
}

// -----------------------------------------------------------------------
void B4bEventAction::getCellSize(){
  auto runData
    = static_cast<B4bRunData*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset(); 

  //   World volume...
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox=dynamic_cast<G4Box*>(worldLV->GetSolid());
  worldZhalf=worldBox->GetZHalfLength(); 

  // Clorimeter volume...
  auto calorLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Calorimeter");
  G4Box* calorBox=dynamic_cast<G4Box*>(calorLV->GetSolid());
  caloXhalf=calorBox->GetXHalfLength();
  caloYhalf=calorBox->GetYHalfLength();
  caloZhalf=calorBox->GetZHalfLength();

  // Layer volume...
  auto layerLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Layer");
  G4Box* layerBox=dynamic_cast<G4Box*>(layerLV->GetSolid());
  caloDz=layerBox->GetZHalfLength()*2.0;  // fult layer thickness, i.e. cell size
//  caloDx=caloDz/4.; //  a cell is cube. 
//  caloDy=caloDz/4.;  
  caloDx=caloDz; //  a cell is cube. 
  caloDy=caloDz;  

  caloNx=((caloXhalf*2.0)/caloDx)+0.01;  // 0.01 to avoid round-off issue
  caloNy=((caloYhalf*2.0)/caloDy)+0.01;  // 0.01 to avoid round-off issue
  caloNz=((caloZhalf*2.0)/caloDz)+0.01;  // 0.01 to avoid round-off issue

  /* 
  std::cout<<"worldZhalf= "<<worldZhalf<<std::endl;
  std::cout<<"caloXhalf = "<<caloXhalf<<"  caloDx "<<caloDx<<"  CaloNx "<<caloNx<<std::endl; 
  std::cout<<"caloYhalf = "<<caloYhalf<<"  caloDy "<<caloDy<<"  CaloNy "<<caloNy<<std::endl; 
  std::cout<<"caloZhalf = "<<caloZhalf<<"  caloDz "<<caloDz<<"  CaloNz "<<caloNz<<std::endl; 
  */

  // save those in TTree...
  
  mNxCell=caloNx;
  mNyCell=caloNy;
  mNzCell=caloNz;
  mDxCell=caloDx;
  mDyCell=caloDy;
  mDzCell=caloDz;
  mzhalfWorld=worldZhalf;

}
