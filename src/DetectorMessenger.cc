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
/// \file electromagnetic/TestEm5/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
  :G4UImessenger(), fTestemDir(0), fDetDir(0), fWinMaterCmd(0), fGasMaterCmd(0), fWinThickCmd(0),
   fGasThickCmd(0), fGasSizYZCmd(0), fGasXposCmd(0), fWorldMaterCmd(0),
   fWorldXCmd(0), fWorldYZCmd(0), fUpdateCmd(0),
   fIonCmd(0), fDetector(Det), fAnodePosCmd(0), fAnodeLengthCmd(0),
   fSegmentLengthCmd(0)
{
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("UI commands specific to this example.");

  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");

  fGasMaterCmd = new G4UIcmdWithAString("/testem/det/setGasMat",this);
  fGasMaterCmd->SetGuidance("Select Material of the Gas.");
  fGasMaterCmd->SetParameterName("choice",false);
  fGasMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWorldMaterCmd = new G4UIcmdWithAString("/testem/det/setWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select Material of the World.");
  fWorldMaterCmd->SetParameterName("wchoice",false);
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWinMaterCmd = new G4UIcmdWithAString("/testem/det/setWinMat", this);
  fWinMaterCmd->SetGuidance("Select Material of the Window.");
  fWinMaterCmd->SetParameterName("wchoice",false);
  fWinMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWinThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWinThick",this);
  fWinThickCmd->SetGuidance("Set Thickness of the Window");
  fWinThickCmd->SetParameterName("SizeZ",false);
  fWinThickCmd->SetRange("SizeZ>0.");
  fWinThickCmd->SetUnitCategory("Length");
  fWinThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGasThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGasThick",this);
  fGasThickCmd->SetGuidance("Set Thickness of the Gas");
  fGasThickCmd->SetParameterName("SizeZ",false);
  fGasThickCmd->SetRange("SizeZ>0.");
  fGasThickCmd->SetUnitCategory("Length");
  fGasThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGasSizYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGasYZ",this);
  fGasSizYZCmd->SetGuidance("Set sizeYZ of the Gas");
  fGasSizYZCmd->SetParameterName("SizeYZ",false);
  fGasSizYZCmd->SetRange("SizeYZ>0.");
  fGasSizYZCmd->SetUnitCategory("Length");
  fGasSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGasXposCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGasXpos",this);
  fGasXposCmd->SetGuidance("Set X pos. of the Gas");
  fGasXposCmd->SetParameterName("Xpos",false);
  fGasXposCmd->SetUnitCategory("Length");
  fGasXposCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWorldXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldX",this);
  fWorldXCmd->SetGuidance("Set X size of the World");
  fWorldXCmd->SetParameterName("WSizeX",false);
  fWorldXCmd->SetRange("WSizeX>0.");
  fWorldXCmd->SetUnitCategory("Length");
  fWorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWorldYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldYZ",this);
  fWorldYZCmd->SetGuidance("Set sizeYZ of the World");
  fWorldYZCmd->SetParameterName("WSizeYZ",false);
  fWorldYZCmd->SetRange("WSizeYZ>0.");
  fWorldYZCmd->SetUnitCategory("Length");
  fWorldYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

  fIonCmd = new G4UIcmdWithADoubleAndUnit("/testem/setPairEnergy", this);
  fIonCmd->SetGuidance("Set energy per electron-ion pair for detector");
  fIonCmd->SetParameterName("en",false,false);
  fIonCmd->SetUnitCategory("Energy");
  fIonCmd->SetDefaultUnit("MeV");
  fIonCmd->SetRange("en>0.");
  fIonCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAnodePosCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAnodePosition", this);
  fAnodePosCmd->SetGuidance("Set anode position with respect to the gas volume");
  fAnodePosCmd->SetParameterName("fAnodePosition", false);
  fAnodePosCmd->SetUnitCategory("Length");
  fAnodePosCmd->SetDefaultUnit("cm");

  fAnodeLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAnodeLength", this);
  fAnodeLengthCmd->SetGuidance("Set the length of the anode volume.");
  fAnodeLengthCmd->SetParameterName("fAnodeX", false);
  fAnodeLengthCmd->SetUnitCategory("Length");
  fAnodeLengthCmd->SetDefaultUnit("cm");
  fAnodeLengthCmd->SetRange("fAnodeX>0.");

  fSegmentLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSegmentLength", this);
  fSegmentLengthCmd->SetGuidance("Set length of each segment within anode volume.");
  fSegmentLengthCmd->SetParameterName("fSegmentX", false);
  fSegmentLengthCmd->SetUnitCategory("Length");
  fSegmentLengthCmd->SetDefaultUnit("cm");
  fSegmentLengthCmd->SetRange("fSegmentX>0.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWinMaterCmd;
  delete fWinThickCmd;
  delete fGasMaterCmd;
  delete fGasThickCmd;
  delete fGasSizYZCmd;
  delete fGasXposCmd;
  delete fWorldMaterCmd;
  delete fWorldXCmd;
  delete fWorldYZCmd;
  delete fUpdateCmd;
  delete fDetDir;
  delete fTestemDir;
  delete fIonCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fWinMaterCmd)
    {fDetector->SetWindowMaterial(newValue);}

  if ( command == fGasMaterCmd )
   {fDetector->SetGasMaterial(newValue);}

  if ( command == fWorldMaterCmd )
   {fDetector->SetWorldMaterial(newValue);}

  if ( command == fGasThickCmd )
  {fDetector->SetGasThickness(fGasThickCmd->GetNewDoubleValue(newValue));}

  if ( command == fWinThickCmd )
  {fDetector->SetWindowThickness(fWinThickCmd->GetNewDoubleValue(newValue));}

  if ( command == fGasSizYZCmd )
   {fDetector->SetGasSizeYZ(fGasSizYZCmd->GetNewDoubleValue(newValue));}

  if ( command == fGasXposCmd )
   {fDetector->SetGasXpos(fGasXposCmd->GetNewDoubleValue(newValue));}

  if ( command == fWorldXCmd )
   {fDetector->SetWorldSizeX(fWorldXCmd->GetNewDoubleValue(newValue));}

  if ( command == fWorldYZCmd )
   {fDetector->SetWorldSizeYZ(fWorldYZCmd->GetNewDoubleValue(newValue));}

  if ( command == fAnodePosCmd )
    {fDetector->SetAnodePosition(fAnodePosCmd->GetNewDoubleValue(newValue));}

  if ( command == fAnodeLengthCmd )
    {fDetector->SetAnodeLength(fAnodeLengthCmd->GetNewDoubleValue(newValue));}

  if ( command == fSegmentLengthCmd )
    {fDetector->SetSegmentLength(fSegmentLengthCmd->GetNewDoubleValue(newValue));}

  if  ( command == fUpdateCmd )
   {fDetector->UpdateGeometry(); }

  if ( command == fIonCmd )
    {
      fDetector->SetPairEnergy(fIonCmd->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
