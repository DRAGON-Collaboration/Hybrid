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
/// \file electromagnetic/Hybrid/Hybrid.cc
/// \brief Main program of the electromagnetic/Hybrid example
//
// $Id: Hybrid.cc 67235 2013-02-08 16:34:49Z vnivanch $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "G4AutoDelete.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  // PhysicsList* physics = new PhysicsList; //Uses PhysicsList.cc


  //Uses reference list QGSP_INCLXX_LIV, EmLivermore Physics with elastic and inelastic hadron scattering physics
  //does not include em options for fluorescence and augor e- specified in physicslist.cc
  G4int verbose = 1;
  G4PhysListFactory factory;
  G4VModularPhysicsList* physics = factory.GetReferencePhysList("QGSP_INCLXX_LIV");
  physics->SetVerboseLevel(verbose);


  runManager->SetUserInitialization(physics);


  // PrimaryGeneratorAction* gun = new PrimaryGeneratorAction();
  //runManager->SetUserInitialization(new DetectorConstruction(gun));

  // set user action classes
  PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  RunAction* run_action = new RunAction(detector, gen_action);
  runManager->SetUserAction(run_action);
  EventAction* event_action = new EventAction(gen_action, run_action, detector);
  runManager->SetUserAction(event_action);
  runManager->SetUserAction(new TrackingAction(detector, run_action, event_action));
  runManager->SetUserAction(new SteppingAction(detector, gen_action, event_action));

  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  else           //define visualization and UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif

#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      UI->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart();
      delete ui;
#endif

#ifdef G4VIS_USE
      delete visManager;
#endif
    }

  // job termination
  //
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
