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
/// \file B4bSteppingAction.cc
/// \brief Implementation of the B4bSteppingAction class

#include "B4bSteppingAction.hh"
#include "B4bRunData.hh"
#include "B4DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "CaloDataStruc.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bSteppingAction::B4bSteppingAction(
                      B4bEventAction* eventAction)
  : G4UserSteppingAction(),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bSteppingAction::~B4bSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bSteppingAction::UserSteppingAction(const G4Step* step)
{
   G4Track* track = step ->GetTrack();
// Collect energy and track length step by step

  // get volume of the current step
  // auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  G4TrackStatus tkstatus=step->GetTrack()->GetTrackStatus();
  const G4DynamicParticle* dynamicParticle= track ->GetDynamicParticle();
  G4int pdgcode=dynamicParticle->GetPDGcode();
  G4int absPdgCode=abs(pdgcode);
  G4ParticleDefinition* particle = dynamicParticle->GetDefinition();
  G4String particleName= particle ->GetParticleName();
  // G4double kinEnergy=dynamicParticle->GetKineticEnergy();

  // if(tkstatus==fStopAndKill && absPdgCode>100 && track->GetTrackID()==1) {
  if(tkstatus==fStopAndKill && absPdgCode>100) {
     fEventAction->FillSecondaries(step);
  }

  //  energy deposit in cell...
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  // auto depth = touchable->GetHistory()->GetDepth();
  auto thisPhysical = touchable->GetVolume(); // mother
  // auto thisCopyNo = thisPhysical->GetCopyNo();
  auto thisName = thisPhysical->GetName();
  G4ThreeVector posA = step->GetPreStepPoint()->GetPosition();

  //  if(thisName.compare(0,5,"Layer")==0)  {
  if(thisName.compare(0,6,"Sensor")==0)  {  
     // std::cout<<"Stepping Action:  volume "<<thisName<<"  copy no "<<std::endl;
     CaloStepData aHit;
     aHit.x=posA.x();
     aHit.y=posA.y();
     aHit.z=posA.z();
     aHit.pid=pdgcode;
     aHit.trackid=track->GetTrackID();
     aHit.globaltime=track->GetGlobalTime();
     aHit.steplength=track->GetTrackLength();
     aHit.edep=edep;

     fEventAction->AccumulateCaloHits(aHit);
  }

  fEventAction->StepAnalysis(step);

}  // end of B4bSteppingAction::UserSteppingAction.


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
