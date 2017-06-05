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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"

#include "G4IonParametrisedLossModel.hh"
#include "G4ParticleDefinition.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4LossTableManager.hh"
#include "G4ElectronIonPair.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, PrimaryGeneratorAction* kin, EventAction* evt)
: G4UserSteppingAction(), 
  fDetector(det), fPrimary(kin), fEventAction(evt), fElIonPair(0)                                         
{
  fElIonPair = G4LossTableManager::Instance()->ElectronIonPair();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{  
   //get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4Material* TrackingMaterial = aStep -> GetTrack() ->
  GetVolume() -> 
  GetLogicalVolume() ->
  GetMaterial(); 
  G4double IonsAlongStep = fElIonPair->SampleNumberOfIonsAlongStep(aStep);

  //Get step particle definition
  G4ParticleDefinition* stepParticle = aStep-> GetTrack() -> GetDefinition();

  //collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double edepNonIon = aStep->GetNonIonizingEnergyDeposit();
  G4double StepLength = aStep->GetStepLength();

  if (volume == fDetector->GetWindow()) fEventAction->AddWindowNonIon(edepNonIon);
  if (volume == fDetector->GetDSSSD()) fEventAction->AddDSSSDNonIon(edepNonIon);
  if (TrackingMaterial == fDetector->GetGasMaterial()) fEventAction->AddGasNonIon(edepNonIon);
      
      
  if (volume == fDetector->GetWindow()) fEventAction->AddWindow(edep);
  if (volume == fDetector->GetDSSSD()) fEventAction->AddDSSSD(edep);
  if (TrackingMaterial == fDetector->GetGasMaterial())
    {
      fEventAction->AddGas(edep);
      fEventAction->AddIons(IonsAlongStep);
      if (stepParticle == fPrimary->GetParticleGun()->GetParticleDefinition())
	{fEventAction->AddStepLength(StepLength);}
      }
 
  //Log Energy Deposition in each anode segment
 
      if (volume->GetCopyNo() != 0) 
	{
	  fEventAction->AddSegment(volume->GetCopyNo(), edep, IonsAlongStep);
				   //, StepLength);
	}
     
 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

