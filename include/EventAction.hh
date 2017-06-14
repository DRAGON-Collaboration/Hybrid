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
/// \file electromagnetic/Si_Ion_Chamber_v7/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 67268 2013-02-13 11:38:40Z ihrivnac $
//
//---------------------------------------------------------------------------
//
// ClassName:   EventAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------

#ifndef EventAction_h
#define EventAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RunAction;
class EventActionMessenger;
class HistoManager;
class DetectorConstruction;

class EventAction : public G4UserEventAction
{
public: // Without description

  EventAction(PrimaryGeneratorAction* kin, RunAction* run, DetectorConstruction* det);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);

    //functions accumulating run variables
  void AddWindow(G4double de) {fEnergyWindow+=de;};
  void AddGas(G4double de) {fEnergyGas += de;};
  void AddDSSSD(G4double de) {fEnergyDSSSD += de;};
  void AddSegment(G4int, G4double, G4double);
          //, G4double);
  void AddIons(G4double de) {totIons += de;};
  void AddStepLength(G4double de) {totStepLength += de;};

  void AddWindowNonIon(G4double de) { fNonIonWindow += de;};
  void AddGasNonIon(G4double de) { fNonIonGas += de;};
  void AddDSSSDNonIon(G4double de) { fNonIonDSSSD += de;};


private:
  PrimaryGeneratorAction* fPrimary;
  G4int fPrintModulo;
  RunAction* fRunAct;
  DetectorConstruction* fDetector;
  HistoManager* fHisto;

  G4int NumberOfSegments;
  G4double fEnergyWindow, fEnergyDSSSD, fEnergyGas;
  G4double fNonIonWindow, fNonIonGas, fNonIonDSSSD;
  G4double IonsAlongStep;
  G4double totIons, totStepLength;

  std::vector<G4double> SegmentEdep;
  std::vector<G4double> SegmentIons;
  std::vector<G4double> SegmentStepLength;


};

#endif
