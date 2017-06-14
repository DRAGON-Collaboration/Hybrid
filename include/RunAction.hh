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
/// \file electromagnetic/TestEm0/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

  void fillPerEvent(G4double, G4double, G4double, std::vector<G4double>&, std::vector<G4double>&);
  void fillPerEventNonIonizing(G4double, G4double, G4double);

  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager* fHisto;

  std::vector<G4double> fSumESegment, fSum2ESegment, fSumSegmentIons, fSum2SegmentIons, Vmax;
  G4double fSumEWindow, fSum2EWindow;
  G4double fSumEDSSSD, fSum2EDSSSD;
  G4double fSumEGas, fSum2EGas;
  G4int NumberOfSegments;

  G4double fSumEWindowNonIon, fSum2EWindowNonIon;
  G4double fSumEDSSSDNonIon, fSum2EDSSSDNonIon;
  G4double fSumEGasNonIon, fSum2EGasNonIon;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
