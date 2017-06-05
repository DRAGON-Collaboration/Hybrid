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
/// \file electromagnetic/Si_Ion_Chamber_v7/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 79241 2014-02-20 16:03:35Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to make analysis and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4StatDouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class G4Step;
class G4ElectronIonPair;
class HistoMessenger;
class G4string;
class G4ParticleDefinition;
class PrimaryGeneratorAction;

class TFile;
class TTree;
class TH1D;
class TH2D;

//const G4int MaxHisto = 7;
const G4int Max2Histo = 3;

class HistoManager
{
public:
  // With description

  static HistoManager* GetPointer();

 // Without description

  HistoManager();
  ~HistoManager();

  void book();
  void save();

  void BeginOfEvent();
  void EndOfEvent(PrimaryGeneratorAction*);
  void BeginOfRun();
  void EndOfRun(G4ParticleDefinition*);
  void SetDetectorParameters(G4int, G4double, G4double, G4double, G4double, G4double, G4double, G4String);
  void fillPerEvent(G4double, G4double, G4double, std::vector<G4double>&, std::vector<G4double>&, G4double, G4double);
  void SetPosition(G4double, G4double, G4double);
  void AddPositionOverflow(){PositionOverflow+=1;};
  void PositionCalc();

  void SetRootFile(G4String);
  G4double CalcStoppingPower(std::vector<G4double>&, G4double);
  void CalcSignal();
  G4double Covariance(std::vector<G4double>&, std::vector<G4double>&);
  G4double Covariance(std::vector<G4double>&, std::vector<G4double>&, G4int);
  G4double Mean(std::vector<G4double>&);
  G4double Mean(std::vector<G4double>&, G4int);
  G4double StdDev(std::vector<G4double>&);
  G4double StdDev(std::vector<G4double>&, G4int);
  std::vector<G4double> PCAa(G4double [2][2], std::vector<G4double>&, std::vector<G4double>&);
  std::vector<G4double> PCAb(G4double [2][2], std::vector<G4double>&, std::vector<G4double>&);
  void CalcInverseP(G4double [2][2]);
  void OutputToTxt(std::vector<G4double>&, std::vector<G4double>&, G4double, G4double, std::vector<G4double>&, std::vector<G4double>&, std::vector<std::vector<G4double> >&, std::vector<std::vector<G4double> >&);
  std::vector<G4double> DiffPulseHeight(std::vector<G4double>&);

  void Fill2Histo(G4int id, G4double ybin, G4double zbin, G4double weight = 1.0);

private:
  G4String fileName;
  static HistoManager* fManager;
  TFile* rootFile;
  TTree* ntupl;
  HistoMessenger* fMessenger;
  G4int EventCounter;
  TH2D* histo2[Max2Histo];
  G4int HeaderFlag, RunHeaderFlag;

  //Detector Parameters
  G4int nSegments;
  G4double SegWidth, GasDensity, IonPairEnergy;
  G4double DSSSDDetectorSizeYZ;
  G4double ChamberLength;
  G4double GasSizeYZ;
  G4double EventEnergy;
  G4String GasName;

  //Energy Deposition as recorded per event
  G4double EdepWindow, EdepDSSSD, EdepGas, IonPair;
  std::vector<G4double> EdepSegment, SegIons, SegStepLength;
  std::vector<G4double> IonPairVector;
  std::vector<G4double> totStepLengthVector;

  //EnergyDeposition as recorded per run
  std::vector<G4double> EdepWindowVector, EdepDSSSDVector, EdepGasVector;
  std::vector<std::vector<G4double> > EdepSegmentVector;
  std::vector<std::vector<G4double> > SegIonsVector;

  //End of track positions
  G4double PositionX, PositionY, PositionZ;
  std::vector<G4double> PositionXVector, PositionYVector, PositionZVector;
  G4int PositionOverflow;
  G4double EventRange;

  //Values calculated at end of run
  G4double InverseP[2][2];
  std::vector<G4double> Vmax, VmaxSTDDEV;
  G4double Efficiency;
  //std::vector<G4double> a, b;
  //std::vector<std::vector<G4double> > aseg, bseg;
};

#endif
