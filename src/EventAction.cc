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
/// \file electromagnetic/Si_Ion_Chamber_v7/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//---------------------------------------------------------------------------
//
// ClassName:   EventAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction(PrimaryGeneratorAction* kin, RunAction* run, DetectorConstruction* det)
  : G4UserEventAction(), fPrimary(kin), fPrintModulo(0), fRunAct(run), fDetector(det), fHisto(HistoManager::GetPointer())
{
	//Print eventID every 'fPrintModulo' event.
    fPrintModulo = 100;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event* evt)
{
	//Clear event data vectors
    SegmentEdep.clear(); 
    SegmentIons.clear();
    SegmentStepLength.clear();
    totIons = totStepLength = 0.;

    NumberOfSegments = fDetector -> GetNumberOfSegments();
    
    //Define size of Vectors containing energy information in each segment
    for (int i=0; i< NumberOfSegments; i++)
      {
	SegmentEdep.push_back(0.);
	SegmentIons.push_back(0.);
	SegmentStepLength.push_back(0.);
      }

  G4int evtNb = evt->GetEventID();
  if (evtNb%fPrintModulo == 0)
    G4cout << "\n--> Begin of event:  " << evtNb << G4endl;

	//Initialize data members
  fEnergyWindow = 0.;
  fEnergyDSSSD = 0.;
  fEnergyGas = 0.;
  fNonIonWindow = 0.;
  fNonIonDSSSD = 0.;
  fNonIonGas =0.;
  fHisto->BeginOfEvent();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::AddSegment(G4int CopyNo, G4double edep, G4double Ions)			     
{
  G4int n = CopyNo-1;
  SegmentEdep[n] += edep;
  SegmentIons[n] += Ions;
}



//...ooo00000ooo........ooo00000ooo.........ooo00000ooo........ooo00000ooo....
void EventAction::EndOfEventAction(const G4Event*)
{
  //accumulates statistic
  fRunAct -> fillPerEvent(fEnergyWindow, fEnergyDSSSD, fEnergyGas, SegmentEdep, SegmentIons);
  fRunAct->fillPerEventNonIonizing(fNonIonWindow, fNonIonDSSSD, fNonIonGas);
  fHisto-> fillPerEvent(fEnergyWindow, fEnergyDSSSD, fEnergyGas, SegmentEdep, SegmentIons, 
			totIons, totStepLength);
  //Histomanager called at end of each event
  fHisto->EndOfEvent(fPrimary);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
