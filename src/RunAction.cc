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
/// \file electromagnetic/TestEm0/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 68243 2013-03-19 18:03:11Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "HistoManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"

#include "G4Run.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Electron.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4IonParametrisedLossModel.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Step.hh"

#include <vector>
#include <iostream>

#include "G4IonParametrisedLossModel.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
  :G4UserRunAction(),fDetector(det), fPrimary(kin), fHisto(HistoManager::GetPointer())
{
  NumberOfSegments = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "RunAction: Begin of run actions are started" << G4endl;
  
  fDetector->PrintCalorParameters();
  
  //Tell HistoManager the parameters of the detector needed to calculate stopping power
  fHisto->SetDetectorParameters(fDetector->GetNumberOfSegments(),
				fDetector->GetSegmentThickness(), 
				fDetector->GetGasMaterial()->GetDensity(),
				fDetector->GetGasMaterial()->GetIonisation()
				->GetMeanEnergyPerIonPair(), 
				fDetector->GetDSSSDDetectorSizeYZ(),
				fDetector->GetGasThickness(),
				fDetector->GetGasSizeYZ(),
				fDetector->GetGasMaterial()->GetName()); 
  //Must come before fHisto->BeginOfRun()
  
  G4int id = aRun->GetRunID();
  G4cout << "### Run " << id << " start" << G4endl;
  
  
  fHisto -> BeginOfRun();
  
  
  //initialize cumulative quantities
  //
  
  fSumESegment.clear();
  fSum2ESegment.clear();
  fSumSegmentIons.clear();
  fSum2SegmentIons.clear();
  Vmax.clear();
  
  NumberOfSegments = fDetector->GetNumberOfSegments();
  
  //Segment data storage vector sizes are defined by the number of segments
  for (int i=0; i < NumberOfSegments; i++)
    {
      fSumESegment.push_back(0.);
      fSum2ESegment.push_back(0.);
      fSumSegmentIons.push_back(0.);
      fSum2SegmentIons.push_back(0.);
      Vmax.push_back(0.);
    }
  
  fSumEWindow = fSum2EWindow = fSumEDSSSD = fSum2EDSSSD = fSumEGas = fSum2EGas = 0.;
  fSumEWindowNonIon = fSum2EWindowNonIon = fSumEDSSSDNonIon = fSum2EDSSSDNonIon = fSumEGasNonIon = fSum2EGasNonIon = 0.;
  
  
  //histograms/ntuples
  //
  fHisto->book(); 
  
#ifdef G4VIS_USE
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  if(pVVisManager)
    {UI->ApplyCommand("/vis/scene/notifyHandlers");}
#endif   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "RunAction: End of run actions are started" << G4endl;
  
  
  // get material
  G4Material* material = fDetector->GetMaterial();
  G4String matName     = material->GetName();
  G4double density     = material->GetDensity();
  G4double radl        = material->GetRadlen();
  
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms for total energy deposition
  //
  fSumEWindow /= NbOfEvents; fSum2EWindow /= NbOfEvents;
  G4double rmsEWindow = fSum2EWindow - fSumEWindow*fSumEWindow;
  if (rmsEWindow >0.) rmsEWindow = std::sqrt(rmsEWindow); else rmsEWindow = 0.;
  
  fSumEDSSSD /= NbOfEvents; fSum2EDSSSD /= NbOfEvents;
  G4double rmsEDSSSD = fSum2EDSSSD - fSumEDSSSD*fSumEDSSSD;
  if (rmsEDSSSD >0.) rmsEDSSSD = std::sqrt(rmsEDSSSD); else rmsEDSSSD = 0.;
  
  fSumEGas /= NbOfEvents; fSum2EGas /= NbOfEvents;
  G4double rmsEGas = fSum2EGas - fSumEGas*fSumEGas;
  if (rmsEGas > 0.) rmsEGas = std::sqrt(rmsEGas); else rmsEGas = 0.;
  
  
  std::vector<G4double> rmsESegment;
  rmsESegment.clear();
  for (unsigned i=0; i < fSumESegment.size(); i++)
    { rmsESegment.push_back(0.);}
  
  std::vector<G4double> rmsSegmentIons;
  rmsSegmentIons.clear();
  for (unsigned i=0; i < fSumSegmentIons.size(); i++)
    { rmsSegmentIons.push_back(0.);}
  
  for (int i=0; i < NumberOfSegments; i++)
    {
      fSumESegment[i] /= NbOfEvents; fSum2ESegment[i] /= NbOfEvents;
      rmsESegment[i] = fSum2ESegment[i] - fSumESegment[i]*fSumESegment[i];
      if (rmsESegment[i] > 0.) rmsESegment[i] = std::sqrt(rmsESegment[i]); else rmsESegment[i] = 0.;
      
      fSumSegmentIons[i] /= NbOfEvents; fSum2SegmentIons[i] /= NbOfEvents;
      rmsSegmentIons[i] = fSum2SegmentIons[i] - fSumSegmentIons[i]*fSumSegmentIons[i];
      if (rmsSegmentIons[i] > 0.) rmsSegmentIons[i] = std::sqrt(rmsSegmentIons[i]); else rmsSegmentIons[i] = 0.;
    }
  
  //Compute statistics for non-ionizing deposition
  fSumEWindowNonIon /= NbOfEvents; fSum2EWindowNonIon /= NbOfEvents;
  G4double rmsEWindowNonIon = fSum2EWindowNonIon - fSumEWindowNonIon*fSumEWindowNonIon;
  if (rmsEWindowNonIon >0.) rmsEWindowNonIon = std::sqrt(rmsEWindowNonIon); else rmsEWindowNonIon = 0.;
  
  fSumEDSSSDNonIon /= NbOfEvents; fSum2EDSSSDNonIon /= NbOfEvents;
  G4double rmsEDSSSDNonIon = fSum2EDSSSDNonIon - fSumEDSSSDNonIon*fSumEDSSSDNonIon;
  if (rmsEDSSSDNonIon >0.) rmsEDSSSDNonIon = std::sqrt(rmsEDSSSDNonIon); else rmsEDSSSDNonIon = 0.;
  
  fSumEGasNonIon /= NbOfEvents; fSum2EGasNonIon /= NbOfEvents;
  G4double rmsEGasNonIon = fSum2EGasNonIon - fSumEGasNonIon*fSumEGasNonIon;
  if (rmsEGasNonIon > 0.) rmsEGasNonIon = std::sqrt(rmsEGasNonIon); else rmsEGasNonIon = 0.;
 
  //print
  //
  G4cout
    << "\n--------------------End of Run------------------------------\n"
    << "\n mean Energy in Window : " << G4BestUnit(fSumEWindow,"Energy")
    << " +- "                          << G4BestUnit(rmsEWindow,"Energy")  
    << "\n mean Energy in DSSSD  : " << G4BestUnit(fSumEDSSSD,"Energy")
    << " +- "                          << G4BestUnit(rmsEDSSSD,"Energy")
    << "\n mean Energy in Gas    : " << G4BestUnit(fSumEGas, "Energy")
    << " +- "                        << G4BestUnit(rmsEGas, "Energy")
    << G4endl;
  
  G4cout
    << "\nNon-ionizing Energy\n"
    << "\n mean non-ionizing Energy in Window : " << G4BestUnit(fSumEWindowNonIon,"Energy")
    << " +- "                          << G4BestUnit(rmsEWindowNonIon,"Energy")  
    << "\n mean non-ionizing Energy in DSSSD  : " << G4BestUnit(fSumEDSSSDNonIon,"Energy")
    << " +- "                          << G4BestUnit(rmsEDSSSDNonIon,"Energy")
    << "\n mean non-ionizing Energy in Gas    : " << G4BestUnit(fSumEGasNonIon, "Energy")
    << " +- "                        << G4BestUnit(rmsEGasNonIon, "Energy")
    << G4endl;
  
  for (int i=0; i < NumberOfSegments; i++)
    {
      G4cout
	<< "\n mean Energy in Segment " << i+1 << " : " << G4BestUnit(fSumESegment[i], "Energy")
	<< " +- "                       << G4BestUnit(rmsESegment[i], "Energy")
	<< "\t Ion Pairs Generated in Segment " << i+1 << " : " << fSumSegmentIons[i]
	<< " +- "                       << rmsSegmentIons[i]
	<< G4endl;
    }
  
  //set precision for printing
  G4int prec = G4cout.precision(6);
  
  // get particle 
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
    ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double charge   = particle->GetPDGCharge();    
  G4double energy   = fPrimary->GetParticleGun()
    ->GetParticleEnergy();
  
  G4cout << "\n " << partName << " ("
         << G4BestUnit(energy,"Energy") << ") in " 
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ";   radiation length: "
         << G4BestUnit(radl,   "Length")       << ")" 
	 <<"Charge " << charge << G4endl;
  
  //save histograms
  fHisto->save();
  //fHisto->EndOfRun(particle, fSumEWindow, rmsEWindow, fSumEGas, rmsEGas, fSumEDSSSD, rmsEDSSSD);
  fHisto->EndOfRun(particle);
  
  // reset default precision
  G4cout.precision(prec);
  
#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
#endif
  
  //Ion Parametrised Loss Model to semi-empirically predict stopping powers
  //G4IonParametrisedLossModel* param = new G4IonParametrisedLossModel();
  //param->PrintDEDXTable(particle, material, 1.97*MeV, 1.97*MeV, 0, false);
  //delete param;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EWin, G4double EDet, G4double EGas, std::vector<G4double>& fSegmentEdep, std::vector<G4double>& fSegmentIons)
{
  //accumulate statistic
  //
  fSumEWindow += EWin; fSum2EWindow += EWin*EWin;
  fSumEDSSSD += EDet; fSum2EDSSSD += EDet*EDet;
  fSumEGas += EGas; fSum2EGas += EGas*EGas;
  
  for (int i=0; i < NumberOfSegments; i++)
    {
      fSumESegment[i] += fSegmentEdep[i]; fSum2ESegment[i] += fSegmentEdep[i] * fSegmentEdep[i];
      fSumSegmentIons[i] += fSegmentIons[i]; fSum2SegmentIons[i] += fSegmentIons[i] * fSegmentIons[i];
    }
}

//.....00000000000........00000000..........00000000...........000000000........

void RunAction::fillPerEventNonIonizing(G4double EWinNonIon, G4double EDetNonIon, G4double EGasNonIon)
{
  //accumulate statistics for non-ionizing energy deposition
  //
  fSumEWindowNonIon += EWinNonIon; fSum2EWindowNonIon += EWinNonIon*EWinNonIon;
  fSumEDSSSDNonIon += EDetNonIon; fSum2EDSSSDNonIon += EDetNonIon*EDetNonIon;
  fSumEGasNonIon += EGasNonIon; fSum2EGasNonIon += EGasNonIon*EGasNonIon;
}
