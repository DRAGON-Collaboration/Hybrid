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
/// \file electromagnetic/Si_Ion_Chamber_v7/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 79241 2014-02-20 16:03:35Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcommand.hh"
#include "G4String.hh"
#include "HistoMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <TTree.h>

#include <iostream>
#include <math.h>

//Histomanager handles not only histrograms but most of the valuable output
//of this application as well as statistical calculations.
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(!fManager) {
    fManager = new HistoManager();
  }
  return fManager;
}

HistoManager::HistoManager() : file0(0)
{
  fileName = "Si_Ion_v8.root";

  fMessenger = new HistoMessenger(this);

  // fMessenger = new HistoMessenger(this);
  RunHeaderFlag = 0;
  //t1
  t1 = 0;

  EdepWindow = 0.;
  EdepDSSSD = 0.;
  EdepGas = 0.;

  //histograms
  for (G4int k=0; k<Max2Histo; k++) histo2[k] = 0;

  //Initialize Detector Parameters
  nSegments = 0;
  SegWidth = 0.;
  GasDensity = 0.;
  IonPairEnergy = 0.;
  GasSizeYZ = 0.;
  DSSSDDetectorSizeYZ = 0.;
  ChamberLength = 0.;
  GasName = "NameNotSet";


}



HistoManager::~HistoManager()
{
  delete fMessenger;
  if (file0 ) delete file0;
}

void HistoManager::book()
{
  //Creating a tree container to handle histograms and t1es
  //This tree is associated to an output file

  file0 = new TFile(fileName,"RECREATE");
  if(!file0) {
    G4cout << " HistoManager::book :"
       << " problem creating the ROOT TFile "
       << G4endl;
    return;
  }

  t1 = new TTree("t1", "Event Data");
  t1->Branch("GasName", &GasName, "GasName/D");
  t1->Branch("EventEnergy", &EventEnergy, "EventEnergy/D");
  t1->Branch("EdepWindow" , &EdepWindow, "EdepWindow/D");
  t1->Branch("EdepDSSSD", &EdepDSSSD, "EdepDSSSD/D");
  t1->Branch("EdepGas", &EdepGas, "EdepGas/D");
  t1->Branch("IonPair", &IonPair, "IonPair/D");
  t1->Branch("EventRange", &EventRange, "EventRange/D");
  for (int i = 0; i < nSegments; i++)
    {
      G4String SegmentString = "EdepSegment";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentString += SegmentNumberString;
      t1->Branch(SegmentString, &EdepSegment[i], "EdepSegment/D");
    }
  for (int i = 0; i < nSegments; i++)
    {
      G4String SegmentString = "SegmentIonPairs";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentString += SegmentNumberString;
      t1->Branch(SegmentString, &SegIons[i], "SegIons/D");
    }
  t1->Branch("PositionX", &PositionX, "PositionX/D");
  t1->Branch("PositionY", &PositionY, "PositionY/D");
  t1->Branch("PositionZ", &PositionZ, "PositionZ/D");

  // 2-D Histos
  histo2[0] = new TH2D("h20", "YZ Position of Particle on DSSSD", 100,
               -5.0*CLHEP::cm, 5.0*CLHEP::cm, 100, -5.0*CLHEP::cm, 5.0*CLHEP::cm);
  histo2[1] = new TH2D("h21", "Edep vs E", 1000, 0, 200*CLHEP::MeV, 1000, 0, 100*CLHEP::MeV);

  G4cout << "\n-----> Histogram file is opened in " << fileName << G4endl;
}

void HistoManager::save()
{
  if (file0) {
    file0->Write();       // Writing the histograms to the file
    file0->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
}

void HistoManager::BeginOfRun()
{
  //Reset objects storing run data for end of run calculations
  EdepWindowVector.clear();
  EdepDSSSDVector.clear();
  EdepGasVector.clear();

  EdepSegmentVector.clear();
  EdepSegmentVector.resize(nSegments);
  SegIonsVector.clear();
  SegIonsVector.resize(nSegments);
  IonPairVector.clear();
  totStepLengthVector.clear();

  EdepSegment.clear();
  for (int i = 0; i < nSegments; i++) {EdepSegment.push_back(0.);}
  SegIons.clear();
  for(int i = 0; i < nSegments; i++) {SegIons.push_back(0.);}

  PositionXVector.clear();
  PositionYVector.clear();
  PositionZVector.clear();
  PositionOverflow = 0;

  EventCounter = 0;

  Vmax.clear();
  VmaxSTDDEV.clear();

  Efficiency = 0.;
  HeaderFlag = 0;

  InverseP[0][0] = InverseP[1][1] = 1.;
  InverseP[1][0] = InverseP[0][1] = 0.;
}

void HistoManager::BeginOfEvent()
{
  //Reset event storage objects
  EventCounter += 1;
  EdepWindow = 0.;
  EdepDSSSD = 0.;
  EdepGas = 0.;
  IonPair = 0.;

  EventRange = 0.;

  for (int i=0; i<nSegments; i++) {EdepSegment[i]=0.;}
  for (int i=0; i<nSegments; i++) {SegIons[i]=0.;}
  PositionX = PositionY = PositionZ = 0.;
}

void HistoManager::EndOfEvent(PrimaryGeneratorAction* fPrimary)
{
    //Retrieves primary particle energy
    EventEnergy = fPrimary->GetParticleGun()->GetParticleEnergy()/CLHEP::MeV;

  if (totStepLengthVector[EventCounter-1]) EventRange = totStepLengthVector[EventCounter-1];
  //Increments root t1 with event values
  if (t1) t1->Fill();

  //Records X and Y position of the end of the particle track
  if ( (std::abs(PositionY) <= DSSSDDetectorSizeYZ/2)&&(std::abs(PositionZ) <= DSSSDDetectorSizeYZ/2) )
    {
      Fill2Histo(1, PositionY, PositionZ);
    }
    //If energy is deposited in the active DSSSD material this is recorded
  if (EdepDSSSD != 0.) Fill2Histo(2, EdepDSSSD, EdepGas, 1);

  //Outputs event data to text file
  std::ofstream EventOutput;
  EventOutput.open("events.out", std::ofstream::out | std::ofstream::app);
  if (HeaderFlag == 0)
    {
      HeaderFlag = 1;
    EventOutput << "Gas\t" << "EventEnergy\t" << "EdepWindow\t" << "EdepDSSSD\t" << "EdepGas\t" << "IonPairs\t" << "EventRange\t";

      for (int i = 0; i < nSegments; i++)
    {
      G4String SegmentString = "EdepSegment";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentString += SegmentNumberString;
      EventOutput << SegmentString << "\t";
    }
      for (int i = 0; i < nSegments; i++)
    {
      G4String SegmentString = "SegmentIonPairs";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentString += SegmentNumberString;
      EventOutput << SegmentString << "\t";
    }
      EventOutput << "PositionX\t" << "PositionY\t" << "PositionZ\t" << G4endl;
    }
  EventOutput << GasName << "\t" << EventEnergy << "\t" << EdepWindow << "\t" << EdepDSSSD << "\t" << EdepGas << "\t" << IonPair << "\t" << EventRange << "\t";
  for (int i = 0; i < nSegments; i++)
    {
      EventOutput << EdepSegment[i] << "\t";
    }
  for (int i = 0; i < nSegments; i++)
    {
      EventOutput << SegIons[i] << "\t";
    }
  EventOutput << PositionX << "\t" << PositionY << "\t" << PositionZ << G4endl;
  EventOutput.close();
}


void HistoManager::EndOfRun(G4ParticleDefinition* particle)
{
  //Get Particle Parameters
  G4double AtomicNumber = particle->GetAtomicNumber();
  G4double AtomicMass = particle->GetAtomicMass();

  //End of run calculations
  CalcSignal();

  G4double MeanStepLength = 0.;
  MeanStepLength = Mean(totStepLengthVector);
  G4double StdDevStepLength = 0.;
  StdDevStepLength = StdDev(totStepLengthVector);

  G4double AverageStoppingPower = 0;
  std::vector<G4double> AverageStoppingPowerSeg;
  std::vector<G4double> StoppingPowerSTDEVSeg;

  AverageStoppingPower = CalcStoppingPower(IonPairVector, MeanStepLength);

  for (int i = 0; i < nSegments; i++)
    {
      AverageStoppingPowerSeg.push_back(CalcStoppingPower(SegIonsVector[i], SegWidth));
      StoppingPowerSTDEVSeg.push_back( StdDev(SegIonsVector[i]) / Mean(SegIonsVector[i]) * CalcStoppingPower(SegIonsVector[i], SegWidth) );
    }

  G4double StoppingPowerSTDEV = 0.;
  StoppingPowerSTDEV = std::pow( (StdDev(IonPairVector) / Mean(IonPairVector) ) , 2.)
    + std::pow( (StdDevStepLength / MeanStepLength) , 2.);
  StoppingPowerSTDEV = std::sqrt(StoppingPowerSTDEV);
  StoppingPowerSTDEV = StoppingPowerSTDEV * AverageStoppingPower;

  //Terminal Output

  G4cout << "Stopping power over whole chamber length: " << AverageStoppingPower / GasDensity
     << " +- " << StoppingPowerSTDEV / GasDensity << " MeV*cm^2*mg^-1"
     << G4endl;

  Efficiency = PositionOverflow;
  Efficiency /= EventCounter;
  Efficiency = 1-Efficiency;

  for ( int i=0; i<nSegments; i++)
    {
      G4cout << "Observable Stopping Power in segment " << i+1 << ":\t" << (AverageStoppingPowerSeg[i]) / GasDensity
         << " +- " << (StoppingPowerSTDEVSeg[i]) / GasDensity << " MeV*cm^2*mg^-1"  << G4endl;
    }

  G4cout << "Mean Positions at end of track:\tX: "
     << Mean(PositionXVector) << " +- " << StdDev(PositionXVector) << " mm"
     << "\tY: " << Mean(PositionYVector) << " +- " << StdDev(PositionYVector) << " mm"
     << "\tZ: " << Mean(PositionZVector) << " +- " << StdDev(PositionZVector) << " mm"
     << G4endl;
  G4cout << "Position Overflow Counter:\t " << PositionOverflow << G4endl;
  G4cout << "Efficiency : " << Efficiency << G4endl;

  //Perform Principal Component Analysis on each non-zero pair (E,dE) in each segment and the whole chamber
  //Input Covariance Matrix, E and dE vectors. Output two vectors representing coordinates a and b which are
  //noncovariant. The points (a,b) are in the basis of the eigenvectors of the covariance matrix, transforming
  //a recoil point into this basis and summing the peak separation in quadrature yields the total two-dimensional
  //covariant peak separation.

  std::vector<G4double> a, b;
  std::vector<std::vector<G4double> > aseg, bseg;
  aseg.resize(nSegments); bseg.resize(nSegments);

  G4double CovMatrix[2][2];
  CovMatrix[0][0] = Covariance(EdepDSSSDVector, EdepDSSSDVector,1)/(StdDev(EdepDSSSDVector,1)*StdDev(EdepDSSSDVector,1));
  CovMatrix[0][1] = Covariance(EdepDSSSDVector, EdepGasVector,1)/(StdDev(EdepDSSSDVector,1)*StdDev(EdepGasVector,1));
  CovMatrix[1][0] = Covariance(EdepGasVector, EdepDSSSDVector,1)/(StdDev(EdepGasVector,1)*StdDev(EdepDSSSDVector,1));
  CovMatrix[1][1] = Covariance(EdepGasVector, EdepGasVector,1)/(StdDev(EdepGasVector,1)*StdDev(EdepGasVector,1));

  a = PCAa(CovMatrix, EdepGasVector, EdepDSSSDVector);
  b = PCAb(CovMatrix, EdepGasVector, EdepDSSSDVector);
  for (int s = 0; s < nSegments; s++)
  {
      CovMatrix[0][0] = Covariance(EdepDSSSDVector, EdepDSSSDVector,1);
      CovMatrix[0][1] = Covariance(EdepDSSSDVector, EdepSegmentVector[s],1);
      CovMatrix[1][0] = Covariance(EdepSegmentVector[s], EdepDSSSDVector,1);
      CovMatrix[1][1] = Covariance(EdepSegmentVector[s], EdepSegmentVector[s],1);
      aseg[s] = PCAa(CovMatrix, EdepSegmentVector[s], EdepDSSSDVector);
      bseg[s] = PCAb(CovMatrix, EdepSegmentVector[s], EdepDSSSDVector);
  }

  //Write to .out text file
  OutputToTxt(AverageStoppingPowerSeg, StoppingPowerSTDEVSeg, AtomicNumber, AtomicMass, a, b, aseg, bseg);


  //Calculate differential pulse height spectrum
    std::vector<std::vector<G4double> > DiffPulseHeightVector;
    std::vector<G4double> BinSize;
    BinSize.resize(nSegments);
    DiffPulseHeightVector.resize(nSegments);
    G4int NBins = 100;
    std::vector<G4double> MaxPulse, MinPulse;
    MaxPulse.resize(nSegments);
    MinPulse.resize(nSegments);
    //Bin energy deposition into NBins
    for (unsigned s = 0; s < MaxPulse.size(); s++) {MaxPulse[s] = 0.; MinPulse[s] = EdepSegmentVector[s][0];}
  for (int s = 0; s < nSegments; s++)
    {
        for (unsigned i = 0; i < EdepSegmentVector[s].size(); i++)
        {
            if (EdepSegmentVector[s][i] > MaxPulse[s]) MaxPulse[s] = EdepSegmentVector[s][i];
            if (EdepSegmentVector[s][i] < MinPulse[s]) MinPulse[s] = EdepSegmentVector[s][i];
        }
    BinSize[s] = (MaxPulse[s]-MinPulse[s])/NBins;
    DiffPulseHeightVector[s] = DiffPulseHeight(EdepSegmentVector[s]);
    }

    //Output differential pulse height spectrum to text
        std::ofstream ResOut;
    ResOut.open("Res.out", std::ofstream::out | std::ofstream::app);
    //ResOut << "Mean Counts +- StdDev" << G4endl;
    for (int s = 0; s < nSegments; s++)
    {
        G4String SegmentString = "Segment ";
        G4int n = s+1;
        G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
        SegmentString += SegmentNumberString;
        ResOut << SegmentString << "\t" << "\t";
    }
    ResOut << G4endl;
    //for (int s = 0; s < nSegments; s++) { ResOut << Mean(DiffPulseHeightVector[s]) << "\t" << StdDev(DiffPulseHeightVector[s]) << "\t";}
    //ResOut << G4endl;
    //ResOut << "Resolution" << G4endl;
    //for (int s = 0; s < nSegments; s++) { ResOut << 2.355*StdDev(DiffPulseHeightVector[s]) / Mean(DiffPulseHeightVector[s]) << "\t" << "\t";}
    //ResOut << G4endl << G4endl;
    for (int i = 0; i < NBins; i++)
    {
        for (int s = 0; s < nSegments; s++)
        {
            ResOut << (MinPulse[s] + (i+1)*BinSize[s]) << "\t" << DiffPulseHeightVector[s][i] << "\t";
        }
        ResOut << G4endl;
    }
    ResOut.close();

}

void HistoManager::SetDetectorParameters(G4int NumberOfSegments, G4double SegmentWidth, G4double Density, G4double MeanEnergyPerIonPair,
G4double DSSSDYZ, G4double GasThickness, G4double GasChamberYZ, G4String GasIdentity)
{

  //Global variables for detector parameters set each run from RunAction
  nSegments = NumberOfSegments;
  SegWidth = SegmentWidth/CLHEP::cm;
  GasDensity = Density/CLHEP::mg * CLHEP::cm3;
  IonPairEnergy = MeanEnergyPerIonPair/CLHEP::MeV;
  DSSSDDetectorSizeYZ = DSSSDYZ/CLHEP::mm;
  ChamberLength = GasThickness/CLHEP::mm;
  GasSizeYZ = GasChamberYZ/CLHEP::mm;
  GasName = GasIdentity;
}

void HistoManager::fillPerEvent(G4double fEnergyWindow, G4double fEnergyDSSSD, G4double fEnergyGas, std::vector<G4double>& SegmentEdep, std::vector<G4double>& SegmentIons, G4double fIon, G4double fStepLength)
{
  //Data from each event if stored in a rootfile or kept in a vector for end of run calculations

  EdepWindow = fEnergyWindow/CLHEP::MeV;
  EdepWindowVector.push_back(fEnergyWindow/CLHEP::MeV);

  EdepDSSSD = fEnergyDSSSD/CLHEP::MeV;
  EdepDSSSDVector.push_back(fEnergyDSSSD/CLHEP::MeV);

  EdepGas = fEnergyGas/CLHEP::MeV;
  EdepGasVector.push_back(fEnergyGas/CLHEP::MeV);

  IonPair = fIon;
  IonPairVector.push_back(fIon);

  totStepLengthVector.push_back(fStepLength/CLHEP::cm);

  for (int i=0; i < nSegments; i++) { EdepSegment[i]=SegmentEdep[i]/CLHEP::MeV;}
  for (int i=0; i < nSegments; i++) { EdepSegmentVector[i].push_back(SegmentEdep[i]/CLHEP::MeV);}

  for (int i=0; i < nSegments; i++) { SegIonsVector[i].push_back(SegmentIons[i]);}
  for (int i=0; i < nSegments; i++) { SegIons[i] = SegmentIons[i];}

}

//End of track positions, tracks not ending within the detector are not accounted for and increment the overflow counter
//All tracks not ending in the DSSSD increment the overflow counter
void HistoManager::SetPosition(G4double Xpos, G4double Ypos, G4double Zpos)
{
  PositionX = Xpos/CLHEP::mm;
  PositionY = Ypos/CLHEP::mm;
  PositionZ = Zpos/CLHEP::mm;

  if (std::abs(Xpos) <= ChamberLength/2+1*CLHEP::mm) {PositionXVector.push_back(Xpos/CLHEP::mm);}
  if (std::abs(Ypos) <= DSSSDDetectorSizeYZ/2) {PositionYVector.push_back(Ypos);}
  if (std::abs(Zpos) <= DSSSDDetectorSizeYZ/2) {PositionZVector.push_back(Zpos);}
}











/////////////////////MATH AND OUTPUT////////////////////////////////////////////////



//Allows declaration/renaming of new root file between runs using HistoMessenger
void HistoManager::SetRootFile(G4String file)
{
  if (file0) delete file0;
  fileName = file;
}

//Calculates average stopping power over whole ion range
G4double HistoManager::CalcStoppingPower(std::vector<G4double>& Ions, G4double Length)
{
  G4double StoppingPower = 0;
  StoppingPower = Mean(Ions) * IonPairEnergy / Length;
  return StoppingPower;
}

void HistoManager::CalcSignal()
{
  //Calculate signal induced by ions generated in each segment
  //Signal Amplitude
  G4double A, d, C = 0;
  A = SegWidth / 100. * GasSizeYZ / 1000. ; //Area of anode segment
  d = 0.01*CLHEP::m; //GasSizeYZ / 1000. / 4.; //distance to Frisch grid from anode
  C = 1.0 * 8.85e-12 * A / d; //Assume permittivity of dielectric (low pressure isobutane) is that of vacuum i.e. k = 1
  for (unsigned i = 0; i < SegIonsVector.size(); i++)
    {
      Vmax.push_back(Mean(SegIonsVector[i]) * 1.602e-19 / C);
      VmaxSTDDEV.push_back(StdDev(SegIonsVector[i]) * 1.602e-19 / C);
    }

}
//Calculate covariance of two vectors
G4double HistoManager::Covariance(std::vector<G4double>& X, std::vector<G4double>& Y)
{
  G4double MeanX, MeanY;
  G4double NonZeroCounter = 0.;
  MeanX = MeanY = 0.;
  if (X.size() == Y.size() )
    {
      for (unsigned i = 0; i < X.size(); i++)
    {
      if ( (X[i] !=0.) && (Y[i]!=0.) ) //Ignored events that do not end in DSSSD
        {
          MeanX += X[i];
          NonZeroCounter += 1.;
          MeanY += Y[i];
        }
    }
      MeanX /= NonZeroCounter;
      MeanY /= NonZeroCounter;
    }
  G4double CovXY = 0;
  if (X.size() == Y.size())
    {
      for (unsigned i = 0; i < X.size(); i++)
      {
          if ( (X[i]!=0.) && (Y[i]!=0.) ) {CovXY += (X[i]-MeanX)*(Y[i]-MeanY);}
      }
      if (NonZeroCounter != 1) CovXY /= NonZeroCounter-1;
    }
  return CovXY;
}

//CountFlag specifies to only consider events with tracks ending in the active DSSSD
G4double HistoManager::Covariance(std::vector<G4double>& X, std::vector<G4double>& Y, G4int CountFlag)
{
    (void)CountFlag; //Compiler won't complain about unused variables
  G4double MeanX, MeanY;
  G4double NonZeroCounter = 0.;
  MeanX = MeanY = 0.;
  if ((X.size() == Y.size()) && (X.size() == EdepDSSSDVector.size()))
    {
        for (unsigned i = 0; i < X.size(); i++)
        {
            if ( (X[i] !=0.) && (Y[i]!=0.) && (EdepDSSSDVector[i]!=0.))

            {
            MeanX += X[i];
            NonZeroCounter += 1.;
            MeanY += Y[i];
            }
        }
      MeanX /= NonZeroCounter;
      MeanY /= NonZeroCounter;
    }
  G4double CovXY = 0;
  if ((X.size() == Y.size()) && (X.size() == EdepDSSSDVector.size()))
    {
      for (unsigned i = 0; i < X.size(); i++)
      {
          if ( (X[i]!=0.) && (Y[i]!=0.) && (EdepDSSSDVector[i]!=0.)) {CovXY += (X[i]-MeanX)*(Y[i]-MeanY);}
      }
      if (NonZeroCounter != 1) CovXY /= NonZeroCounter-1;
    }
  return CovXY;
}

G4double HistoManager::Mean(std::vector<G4double>& X)
{
  G4double MeanX = 0.;
  G4double NonZeroCounter = 0.;
  for (unsigned i = 0; i < X.size(); i++)
    {
        if ( (X[i] != 0.) ) {MeanX += X[i]; NonZeroCounter += 1.;}
    }
  if (NonZeroCounter != 0.) MeanX /= NonZeroCounter;
  return MeanX;
}

G4double HistoManager::Mean(std::vector<G4double>& X, G4int CountFlag)
{
    (void)CountFlag;
  G4double MeanX = 0.;
  G4double NonZeroCounter = 0.;
  if (X.size() == EdepDSSSDVector.size())
  {
        for (unsigned i = 0; i < X.size(); i++)
        {
        if ( (X[i] != 0.) && (EdepDSSSDVector[i]!=0.)) {MeanX += X[i]; NonZeroCounter += 1.;}
        }
    }
  if (NonZeroCounter != 0.) MeanX /= NonZeroCounter;
  return MeanX;
}

G4double HistoManager::StdDev(std::vector<G4double>& X)
{
  G4double Variance = 0.;
  G4double NonZeroCounter = 0.;
  for (unsigned i = 0; i < X.size(); i++)
    {
      if (X[i]!=0.) {Variance += std::pow( (X[i]-Mean(X)) , 2.); NonZeroCounter += 1.;}
    }
  if (NonZeroCounter != 1) Variance /= (NonZeroCounter-1.);
  G4double STDDEV;
  STDDEV = std::sqrt(Variance);
  return STDDEV;
}

G4double HistoManager::StdDev(std::vector<G4double>& X, G4int CountFlag)
{
    (void)CountFlag;
  G4double Variance = 0.;
  G4double NonZeroCounter = 0.;
  if (X.size() == EdepDSSSDVector.size())
  {
  for (unsigned i = 0; i < X.size(); i++)
    {
      if ((X[i]!=0) && (EdepDSSSDVector[i]!=0.)) {Variance += std::pow( (X[i]-Mean(X)) , 2.); NonZeroCounter += 1.;}
    }
  }
  if (NonZeroCounter != 1) Variance /= (NonZeroCounter-1.);
  G4double STDDEV;
  STDDEV = std::sqrt(Variance);
  return STDDEV;
}

void HistoManager::Fill2Histo(G4int ih, G4double ybin, G4double zbin, G4double weight)
{
  if (ih >= Max2Histo) {
    G4cout << "---> warning from HistoManager::Fill2Histo() : histo " << ih
           << " does not exist. (ybin=" << ybin << "zbin=" << zbin <<" weight=" << weight << ")"
           << G4endl;
    return;
  }
  if  (histo2[ih]) { histo2[ih]->Fill(ybin, zbin, weight); }

}

std::vector<G4double> HistoManager::PCAa(G4double CovMatrix[2][2], std::vector<G4double>& EGas, std::vector<G4double>& EDSSSD)
{
  CalcInverseP(CovMatrix);
  G4double MeanEdepGas = Mean(EGas,1);
  G4double MeanEdepDSSSD = Mean(EDSSSD,1);
  std::vector<G4double> returnvector;
  if (EGas.size() == EDSSSD.size())
    {
      for (unsigned i=0; i < EGas.size(); i++)
    {
      if ((EGas[i] != 0.) && (EDSSSD[i] != 0.))
        {
          returnvector.push_back( (EDSSSD[i]-MeanEdepDSSSD)*InverseP[0][0] + (EGas[i]-MeanEdepGas)*InverseP[0][1] );
        }
    }
    }
  return returnvector;
}

std::vector<G4double> HistoManager::PCAb(G4double CovMatrix[2][2], std::vector<G4double>& EGas, std::vector<G4double>& EDSSSD)
{
  CalcInverseP(CovMatrix);
  G4double MeanEdepGas = Mean(EGas,1);
  G4double MeanEdepDSSSD = Mean(EDSSSD,1);
  std::vector<G4double> returnvector;
  if (EGas.size() == EDSSSD.size())
    {
      for (unsigned i=0; i < EGas.size(); i++)
    {
      if ((EGas[i] != 0.) && (EDSSSD[i] != 0.))
        {
          returnvector.push_back( (EDSSSD[i]-MeanEdepDSSSD)*InverseP[1][0] + (EGas[i]-MeanEdepGas)*InverseP[1][1] );
        }
    }
    }
  return returnvector;
}

void HistoManager::CalcInverseP(G4double CovMatrix[2][2])
{
  G4double T =  CovMatrix[0][0] + CovMatrix[1][1];
  G4double D = CovMatrix[0][0]*CovMatrix[1][1]-CovMatrix[0][1]*CovMatrix[1][0];
  G4double E1 = T/2. + std::sqrt((std::pow(T,2.)/4-D));
  G4double E2 = T/2. - std::sqrt((std::pow(T,2.)/4-D));
  G4double V1[2][1], V2[2][1];

  if (CovMatrix[0][1] != 0.)
    {
      V1[0][0] = CovMatrix[0][1];
      V1[1][0] = E1 - CovMatrix[0][0];
      V2[0][0] = CovMatrix[0][1];
      V2[1][0] = E2 - CovMatrix[0][0];
    }

  else if (CovMatrix[1][0] != 0.)
    {
      V1[0][0] = E1 - CovMatrix[1][1];
      V1[1][0] = CovMatrix[1][0];
      V2[0][0] = E2 - CovMatrix[1][1];
      V2[1][0] = CovMatrix[1][0];
    }

  else
    {
      V1[0][0] = 1.; V1[1][0] = 0.;
      V2[0][0] = 1.; V2[1][0] = 0.;
    }

  G4double NormV1 = std::sqrt(std::pow(V1[0][0],2.) + std::pow(V1[1][0],2.));
  V1[0][0] = V1[0][0] / NormV1; V1[1][0] = V1[1][0] / NormV1;
  G4double NormV2 = std::sqrt(std::pow(V2[0][0],2.) + std::pow(V2[1][0],2.));
  V2[0][0] = V2[0][0] / NormV2; V2[1][0] = V2[1][0] / NormV2;

  G4double detP;
  detP = V1[0][0]*V2[1][0]-V1[1][0]*V2[0][0];
  InverseP[0][0] = 1./detP * V2[1][0];
  InverseP[0][1] = 1./detP * (-1.) * V2[0][0];
  InverseP[1][0] = 1./detP * (-1.) * V1[1][0];
  InverseP[1][1] = 1./detP * V1[0][0];
}

void HistoManager::OutputToTxt(std::vector<G4double>& AverageStoppingPowerSeg, std::vector<G4double>& StoppingPowerSTDEVSeg, G4double AtomicNumber, G4double AtomicMass,
                   std::vector<G4double>& a, std::vector<G4double>& b, std::vector<std::vector<G4double> >& aseg, std::vector<std::vector<G4double> >& bseg)
{
  (void) a; (void) b; (void) aseg; (void) bseg; (void) AverageStoppingPowerSeg; (void) StoppingPowerSTDEVSeg; //Cast as void to avoid compiler warnings. Comment if you want to use these.
  std::ofstream RunOutput;
  RunOutput.open("run.out", std::ofstream::out | std::ofstream::app);
  //Output Run Headers
  if (RunHeaderFlag == 0)
    {
      RunHeaderFlag = 1;
      RunOutput << "Gas\t" << "Z\t" << "A\t" << "Range (mm)\t" << "Range Straggling (mm)\t";
      RunOutput << "EdepWindow (MeV)\t" << "EdepWindowSigma (MeV)\t" << "EdepGas (MeV)\t" << "EdepGasSigma (MeV)\t" << "EdepDSSSD (MeV)\t" << "EdepDSSSDSigma (MeV)\t";
//		<< "Cov(E,E) (MeV)\t" << "Cov(E,dE) (MeV)\t" << "Cov(dE,E) (MeV)\t" << "Cov(dE,dE) (MeV)\t";
      RunOutput << "IonGas\t" << "IonGasSigma\t" << "Efficiency\t"; //<< "a\t" << "StdDev a\t" << "b\t" << "StdDev b\t";

      for ( int i = 0; i < nSegments; i++)
    {
      G4String SegmentStringA = "EdepSegment (MeV) ";
      G4String SegmentStringB = "EdepSegmentSigma (MeV) ";
//    G4String SegmentStringC = "Cov(E,E) (MeV) ";
//    G4String SegmentStringD = "Cov(E,dE) (MeV) ";
//    G4String SegmentStringE = "Cov(dE,E) (MeV) ";
//    G4String SegmentStringF = "Cov(dE,dE) (MeV) ";
//    G4String SegmentStringG = "Mean a (MeV) ";
//    G4String SegmentStringH = "StdDev a (MeV) ";
//    G4String SegmentStringI = "Mean b (MeV) ";
//    G4String SegmentStringJ = "StdDev b (MeV) ";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentStringA += SegmentNumberString;
      SegmentStringB += SegmentNumberString;
      //SegmentStringC += SegmentNumberString;
      //SegmentStringD += SegmentNumberString;
      //SegmentStringE += SegmentNumberString;
      //SegmentStringF += SegmentNumberString;
//    SegmentStringG += SegmentNumberString;
//    SegmentStringH += SegmentNumberString;
//    SegmentStringI += SegmentNumberString;
//    SegmentStringJ += SegmentNumberString;
      RunOutput << SegmentStringA << "\t" << SegmentStringB << "\t"; //<< SegmentStringC << "\t"
            //<< SegmentStringD << "\t" << SegmentStringE << "\t" << SegmentStringF << "\t"
//          << SegmentStringG << "\t" << SegmentStringH << "\t" << SegmentStringI << "\t"
//          << SegmentStringJ << "\t";
    }
      for ( int i = 0; i < nSegments; i++)
    {
      G4String SegmentStringA = "IonPairSegment";
      G4String SegmentStringB = "IonPairSegmentSigma";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentStringA += SegmentNumberString;
      SegmentStringB += SegmentNumberString;
      RunOutput << SegmentStringA << "\t" << SegmentStringB << "\t";
    }
      /*for ( int i = 0; i < nSegments; i++)
    {
      G4String SegmentStringA = "StoppingPowerSegment (MeVcm2/mg)";
      G4String SegmentStringB = "StoppingPowerSegmentSigma (MeVcm2/mg)";
      G4int n = i+1;
      G4String SegmentNumberString = G4UIcommand::ConvertToString(n);
      SegmentStringA += SegmentNumberString;
      SegmentStringB += SegmentNumberString;
      RunOutput << SegmentStringA << "\t" << SegmentStringB << "\t";
    }*/
      RunOutput << G4endl;
    }

  //Output Run Data
  RunOutput << GasName << "\t" << AtomicNumber << "\t" << AtomicMass << "\t" << Mean(PositionXVector) << "\t" << StdDev(PositionXVector) << "\t";
  RunOutput << Mean(EdepWindowVector) << "\t" << StdDev(EdepWindowVector) << "\t"
        << Mean(EdepGasVector) << "\t" << StdDev(EdepGasVector) << "\t"
        << Mean(EdepDSSSDVector) << "\t" << StdDev(EdepDSSSDVector) << "\t";
       /* << Covariance(EdepDSSSDVector, EdepDSSSDVector) << "\t"
        << Covariance(EdepDSSSDVector, EdepGasVector) << "\t"
        << Covariance(EdepGasVector, EdepDSSSDVector) << "\t"
        << Covariance(EdepGasVector, EdepGasVector) << "\t";*/
  RunOutput << Mean(IonPairVector) << "\t" << StdDev(IonPairVector) << "\t" << Efficiency << "\t";
  //RunOutput << Mean(a) << "\t" << StdDev(a) << "\t" << Mean(b) << "\t" << StdDev(b) << "\t";
  /*for (int i = 0; i < nSegments; i++)
    {
      RunOutput << Vmax[i] << "\t" << VmaxSTDDEV[i] << "\t";
    }*/
  for (int i = 0; i < nSegments; i++)
    {
      RunOutput << Mean(EdepSegmentVector[i]) << "\t" << StdDev(EdepSegmentVector[i]) << "\t";
/*		<< Covariance(EdepDSSSDVector, EdepDSSSDVector) << "\t"
        << Covariance(EdepDSSSDVector, EdepSegmentVector[i]) << "\t"
        << Covariance(EdepSegmentVector[i], EdepDSSSDVector) << "\t"
        << Covariance(EdepSegmentVector[i], EdepSegmentVector[i]) << "\t"
        << Mean(aseg[i]) << "\t" << StdDev(aseg[i]) << "\t"
        << Mean(bseg[i]) << "\t" << StdDev(bseg[i]) << "\t";*/
    }
  for (int i = 0; i < nSegments; i++)
    {
      RunOutput << Mean(SegIonsVector[i]) << "\t" << StdDev(SegIonsVector[i]) << "\t";
    }

  RunOutput << G4endl;
  RunOutput.close();
}

std::vector<G4double> HistoManager::DiffPulseHeight(std::vector<G4double>& pulse)
{
  G4int NBins;
  G4double BinSize;
  G4double MaxPulse, MinPulse;
  MaxPulse = 0.;
  MinPulse = 0.;
  NBins = 100;
  std::vector<G4double> PulseHeight;
  PulseHeight.resize(NBins);
  for (unsigned i = 0; i < pulse.size(); i++)
    {
      if (pulse[i] > MaxPulse) MaxPulse = pulse[i];
      if (i==0) MinPulse = pulse[i];
      if (pulse[i] < MinPulse) MinPulse = pulse[i];
    }
  BinSize = (MaxPulse-MinPulse)/NBins;

  for (int b = 0; b < NBins; b++)
    {
      for (unsigned i = 0; i < pulse.size(); i++)
        {
        if ( ( (MinPulse + b*BinSize) < pulse[i]) && (pulse[i] <= (MinPulse + (b+1)*BinSize)) ) PulseHeight[b] +=1;
        }
    }

  for (unsigned i = 0; i < PulseHeight.size(); i++)
    {
      PulseHeight[i] /= BinSize;
    }
  return PulseHeight;
}
