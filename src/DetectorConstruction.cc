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
/// Si_Ion_Chamber_v8																																																											/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TargetSD.hh"
#include "G4SDManager.hh"
#include "F02ElectricFieldSetup.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4AutoDelete.hh"
//#include "SegParameterisation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fSolidDSSSDActive(0),fLogicDSSSDActive(0),fPhysiDSSSDActive(0),fDSSSDActiveMaterial(0),
 fSolidDSSSDDeadlayer(0),fLogicDSSSDDeadlayer(0),fPhysiDSSSDDeadlayer(0),fDSSSDDeadlayerMaterial(0),
 fSolidDSSSDDetector(0),fLogicDSSSDDetector(0),fPhysiDSSSDDetector(0),fDSSSDDetectorMaterial(0),
 fSolidWindow(0),fLogicWindow(0),fPhysiWindow(0),fWindowMaterial(0),
 fSolidGas(0),fLogicGas(0),fPhysiGas(0),fGasMaterial(0),
  fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),fWorldMaterial(0),fRegGasDet(0),fGasDetectorCuts(0),fTargetSD(0), fDefaultWorld(true)
{
  //Define Detector Geometry
  fDSSSDActiveThickness    = 300*um;
  fDSSSDDetectorSizeYZ     = 49.5*mm;
  fDSSSDActiveSizeZ        = fDSSSDDetectorSizeYZ; // active silicon width
  fDSSSDDeadlayerThickness = 0.4*um;
  fWindowThickness         = 0.5*um;
  fWindowSizeYZ            = 5.0*cm;
  fXposWindow              = 0.0*cm;
  fGasPressure             = 10.0;
  fGasThickness            = 10.8*cm;
  fGasSizeYZ               = 10.0*cm;
  fXposGas                 = 0.0*cm;
  ComputeCalorParameters();

   //Anode Geometry
  fAnodeX = fGasThickness;
  fAnodePosition = 0; //Anode size and position can be set from DetectorMessenger
  fSegmentX = fGasThickness/11; //Segment width can be set from DetectorMessenger

  // materials
  DefineMaterials();
  SetWorldMaterial("Galactic");
  SetGasMaterial("isobutane10torr");
  SetWindowMaterial("G4_MYLAR");
  fDSSSDActiveMaterial    = G4Material::GetMaterial("Silicon");
  fDSSSDDeadlayerMaterial = G4Material::GetMaterial("Aluminium");
  fDSSSDDetectorMaterial  = G4Material::GetMaterial("Galactic");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);

  G4double cut =0.010*mm;
  fGasDetectorCuts = new G4ProductionCuts();
  fGasDetectorCuts->SetProductionCut(cut,"gamma");
  fGasDetectorCuts->SetProductionCut(cut,"e-");
  fGasDetectorCuts->SetProductionCut(cut,"e+");
  fGasDetectorCuts->SetProductionCut(cut,"proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fGasDetectorCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials

  G4String symbol, name;             //a  = molar mass;
  G4double a, z, density;            //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

  //
  // define simple materials
  //

  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr =
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //

  G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* CH = new G4Material("Plastic", density= 1.04*g/cm3, ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  G4Material* PMMA = new G4Material("PMMA", density = 1.190*g/cm3, ncomponents=3);
  PMMA->AddElement(C, natoms = 3);
  PMMA->AddElement(O, natoms = 2);
  PMMA->AddElement(H, natoms = 8);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Air10atm =
    new G4Material("Air10atm", density= 12.93*kg/m3, ncomponents=2,
                   kStateGas, 273.15*kelvin, 10.*atmosphere);
  Air10atm->AddElement(N, fractionmass=0.7);
  Air10atm->AddElement(O, fractionmass=0.3);

  //Graphite
  //
  G4Material* Graphite =
    new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);

  //Havar
  //
  G4Element* Cr = new G4Element("Chrome", "Cr", z=25, a=  51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", z=26, a=  55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", z=27, a=  58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28, a=  58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", z=74, a= 183.850*g/mole);

  G4Material* Havar =
    new G4Material("Havar", density= 8.3*g/cm3, ncomponents=5);
  Havar->AddElement(Cr, fractionmass=0.1785);
  Havar->AddElement(Fe, fractionmass=0.1822);
  Havar->AddElement(Co, fractionmass=0.4452);
  Havar->AddElement(Ni, fractionmass=0.1310);
  Havar->AddElement(W , fractionmass=0.0631);

  //
  // examples of gas
  //
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);

  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);

  G4Material* CO2 =
    new G4Material("CarbonicGas", density= 1.15*mg/cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* ArCO2 =
    new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
  ArCO2->AddElement (Ar,  fractionmass=0.7844);
  ArCO2->AddMaterial(CO2, fractionmass=0.2156);

  new G4Material("NitrogenGas", z=7, a=14.01*g/mole, density=5697.4*kg/m3,
                 kStateGas, 300*kelvin, 10*atmosphere);

  //another way to define mixture of gas per volume
  G4Material* NewArCO2 =
    new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
  NewArCO2->AddElement (Ar, natoms=8);
  NewArCO2->AddElement (C,  natoms=2);
  NewArCO2->AddElement (O,  natoms=4);

  G4Material* ArCH4 =
    new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
  ArCH4->AddElement (Ar, natoms=93);
  ArCH4->AddElement (C,  natoms=7);
  ArCH4->AddElement (H,  natoms=28);

  G4Material* XeCH =
    new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                   kStateGas, 293.15*kelvin, 1*atmosphere);
  XeCH->AddElement (Xe, natoms=875);
  XeCH->AddElement (C,  natoms=225);
  XeCH->AddElement (H,  natoms=700);

  G4Material* steam =
    new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);

  //Isobutane definitions
  G4double IsobutaneTemperature = 273.15*kelvin;
  G4double IsobutaneDensity;
  G4double IsobutanePressure;

  IsobutanePressure = fGasPressure/760.*atmosphere;
  IsobutaneDensity = 0.00341*mg/cm3;

  G4Material* Isobutane = new G4Material(name = "Isobutane", IsobutaneDensity, ncomponents = 2, kStateGas, IsobutaneTemperature, IsobutanePressure);
  Isobutane-> AddElement(C,4);
  Isobutane-> AddElement(H,10);
  Isobutane-> GetIonisation() -> SetMeanExcitationEnergy(0.0000483);
  Isobutane-> GetIonisation() -> SetMeanEnergyPerIonPair(23.0*eV);

  //
  // example of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);

  //Set the Gas Material

  //SetGasMaterial(Isobutane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter

  fXstartGas = fXposGas-0.5*fGasThickness;
  fXendGas   = fXposGas+0.5*fGasThickness;
  fDSSSDDetectorThickness =   fDSSSDActiveThickness + fDSSSDDeadlayerThickness;

  fXposWindow = fXstartGas-0.5*fWindowThickness;
  fXposDSSSDDetector = fXendGas+0.5*fDSSSDDetectorThickness;

  if (fDefaultWorld) {
     fWorldSizeX = 3*m; fWorldSizeYZ= 1.2*fGasSizeYZ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();

  if(fRegGasDet) {delete fRegGasDet; }
  fRegGasDet = new G4Region("GasDetector");
  fRegGasDet->SetProductionCuts(fGasDetectorCuts);

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();

  // World
  //
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);   //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,                //its solid
                                   fWorldMaterial,        //its material
                                   "World");                //its name

  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                                 fLogicWorld,                //its logical volume
                                 "World",                //its name
                                 0,                        //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number

  //Gas
  //
  fSolidGas = new G4Box("Gas",
                      fGasThickness/2,fGasSizeYZ/2,fGasSizeYZ/2);

  fLogicGas = new G4LogicalVolume(fSolidGas,    //its solid
                                            fGasMaterial, //its material
                                          "Gas");       //its name

  fPhysiGas = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(fXposGas,0.,0.),    //its position
                                fLogicGas,     //its logical volume
                                "Gas",         //its name
                                fLogicWorld,        //its mother
                                false,             //no boulean operat
                                0);                //copy number

  //Anode Segments Holder Volume
  //fAnodeX = 30*cm;
  //fAnodePosition = 0; //Anode size and position can be set from DetectorMessenger

  fSolidAnodeBox = new G4Box("Anode Holder Volume", fAnodeX/2, fGasSizeYZ/2, fGasSizeYZ/2);
  fLogicAnodeBox = new G4LogicalVolume(fSolidAnodeBox, fGasMaterial, "Anode Dead Region");
  fPhysiAnodeBox = new G4PVPlacement(0, G4ThreeVector(fAnodePosition,0.,0.),
                     fLogicAnodeBox,
                     "Anode Box",
                     fLogicGas,
                     false,
                     0);
  //Anode Segments
  //fSegmentX = 1.0*cm; //Segment width can be set from DetectorMessenger
  nSegments = floor(fAnodeX / fSegmentX);

  if (nSegments != 1)
    fAnodeDeadLayerThickness = (fAnodeX - (nSegments * fSegmentX)) / (nSegments-1);

  if (nSegments == 1)
    fAnodeDeadLayerThickness = 0*cm;

  fSolidSegment = new G4Box("GasSegment",fSegmentX/2, fGasSizeYZ/2, fGasSizeYZ/2);
  fLogicSegment = new G4LogicalVolume(fSolidSegment, fGasMaterial, "AnodeSegment");

  for (int i = 0; i < nSegments; i++)
    {
      // G4int n = i+1;
      // G4String SegmentNumber = G4UIcommand::ConvertToString(n);

      //      fPhysiSegment =
new G4PVPlacement(0, G4ThreeVector(-fAnodeX/2+fSegmentX/2+i*(fSegmentX+fAnodeDeadLayerThickness), 0., 0.),
                    fLogicSegment,
                    "Gas Segment",
                    fLogicAnodeBox,
                    false,
                    i+1);
                    }

  fRegGasDet->AddRootLogicalVolume(fLogicSegment);
  fRegGasDet->AddRootLogicalVolume(fLogicGas);


   // Ionization Chamber Entrance Window
  //
  fSolidWindow = new G4Box("Window",
                      fWindowThickness/2,fWindowSizeYZ/2,fWindowSizeYZ/2);

  fLogicWindow = new G4LogicalVolume(fSolidWindow,    //its solid
                                            fWindowMaterial, //its material
                                          "Window");       //its name

  fPhysiWindow = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(fXposWindow,0.,0.),    //its position
                                fLogicWindow,     //its logical volume
                                "Window",         //its name
                                fLogicWorld,        //its mother
                                false,             //no boolean operation
                                0);                //copy number

  // DSSSD holder volume for strips and deadlayer
  //
  fSolidDSSSDDetector = new G4Box("DSSSDDetector",
                      fDSSSDDetectorThickness/2,fDSSSDDetectorSizeYZ/2,fDSSSDActiveSizeZ/2);

  fLogicDSSSDDetector = new G4LogicalVolume(fSolidDSSSDDetector,    //its solid
                                          fDSSSDDetectorMaterial, //its material
                                          "DSSSDDetector");       //its name

  fPhysiDSSSDDetector = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(fXposDSSSDDetector,0.,0.),         //its position
                                fLogicDSSSDDetector,             //its logical volume
                                "DSSSDDetector",                 //its name
                                fLogicWorld,         //its mother
                                false,                       //no boolean operator
                                0);                      //copy number
  // DSSSD Active Strip
  //
  fSolidDSSSDActive = new G4Box("DSSSDActive",
                      fDSSSDActiveThickness/2,fDSSSDDetectorSizeYZ/2,fDSSSDActiveSizeZ/2);

  fLogicDSSSDActive = new G4LogicalVolume(fSolidDSSSDActive,    //its solid
                                          fDSSSDActiveMaterial, //its material
                                          "DSSSDActive");       //its name

  fPhysiDSSSDActive = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(fDSSSDDeadlayerThickness/2,0.,0.),         //its position
                                fLogicDSSSDActive,           //its logical volume
                                "DSSSDActive",               //its name
                                fLogicDSSSDDetector,         //its mother
                                false,                       //no boolean operator
                                0);                      //copy number

  // DSSSD Deadlayer
  //
  fSolidDSSSDDeadlayer = new G4Box("DSSSDDeadlayer",
                      fDSSSDDeadlayerThickness/2,fDSSSDDetectorSizeYZ/2,fDSSSDDetectorSizeYZ/2);

  fLogicDSSSDDeadlayer = new G4LogicalVolume(fSolidDSSSDDeadlayer,    //its solid
                                          fDSSSDDeadlayerMaterial, //its material
                                          "DSSSDDeadlayer");       //its name

  fPhysiDSSSDDeadlayer = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(-fDSSSDActiveThickness/2,0.,0.),       //its position
                                fLogicDSSSDDeadlayer,            //its logical volume
                                "DSSSDDeadlayer",                //its name
                                fLogicDSSSDDetector,                 //its mother
                                false,                       //no boolean operator
                                0);                      //copy number


  //Sensitive Detectors
  G4SDManager* SDman =  G4SDManager::GetSDMpointer();
  if(!fTargetSD)
    {
      fTargetSD = new TargetSD("GasSD");
      SDman->AddNewDetector(fTargetSD);
    }
  fLogicSegment->SetSensitiveDetector(fTargetSD);
  fLogicGas->SetSensitiveDetector(fTargetSD);

//Construct the field creator - this will register the field it creates
if (!fEmFieldSetup.Get()) {
    F02ElectricFieldSetup* fieldSetup = new F02ElectricFieldSetup();
    G4AutoDelete::Register(fieldSetup); //Kernel will delete the messenger
    fEmFieldSetup.Put(fieldSetup);}

 PrintCalorParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void DetectorConstruction::ConstructField()
{
//Construct the field creator - this will register the field it creates
if (!fEmFieldSetup.Get()) {
    F02ElectricFieldSetup* fieldSetup = new F02ElectricFieldSetup();
    G4AutoDelete::Register(fieldSetup); //Kernel will delete the messenger
    fEmFieldSetup.Put(fieldSetup);}
    }*/


void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << fGasMaterial << G4endl;
  G4cout << "\n" << fWindowMaterial << G4endl;
  G4cout << "\n" << fDSSSDActiveMaterial << G4endl;

  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeX,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is "
         << G4BestUnit(fWorldSizeYZ,"Length") << G4endl;
  G4cout << "\n The gas Gas is made of "
         <<G4BestUnit(fGasThickness,"Length")
         << " of " << fGasMaterial->GetName();
  G4cout << ". The transverse size (YZ) is "
         << G4BestUnit(fGasSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the Gas "
         << G4BestUnit(fXposGas,"Length");
  G4cout << "\n The anode is positioned at " << G4BestUnit(fAnodePosition, "Length")
     << " with respect to the middle of the gas.";
  G4cout << "\n The anode is " << G4BestUnit(fAnodeX, "Length")
     << " long and divided into " << nSegments << " segments." << G4endl;
  G4cout << "Each segment has a width of " << G4BestUnit(fSegmentX, "Length") << G4endl;

  G4cout << "\n The WINDOW is made of "
         <<G4BestUnit(fWindowThickness,"Length")
         << " of " << fWindowMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the window is "
         << G4BestUnit(fWindowSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the window is "
         << G4BestUnit(fXposWindow,"Length");
  G4cout << "\n The ACTIVE STRIP is made of "
         <<G4BestUnit(fDSSSDActiveThickness,"Length")
         << " of " << fDSSSDActiveMaterial->GetName();
  G4cout << ". The Length of the strip is "
         << G4BestUnit(fDSSSDDetectorSizeYZ,"Length") << G4endl;
  G4cout << ". The Width of the strip is "
         << G4BestUnit(fDSSSDActiveSizeZ,"Length") << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWindowMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWindowMaterial != pttoMaterial) {
    fWindowMaterial = pttoMaterial;
    if(fLogicWindow) fLogicWindow->SetMaterial(fWindowMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetGasMaterial(
                      //G4Material* materialChoice
G4String materialChoice
)
{
  // search the material by its name
  G4Material* pttoMaterial =
    // materialChoice;
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fGasMaterial != pttoMaterial) {
    fGasMaterial = pttoMaterial;
    if(fLogicGas) fLogicGas->SetMaterial(fGasMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWindowThickness(G4double val)
{
  fWindowThickness = val;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasPressure(G4double val)
{
  fGasPressure = val;
  fAnodeX = fGasThickness;
  fSegmentX = fGasThickness/10;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasThickness(G4double val)
{
  fGasThickness = val;
  fAnodeX = fGasThickness;
  fSegmentX = fGasThickness/10;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasSizeYZ(G4double val)
{
  fGasSizeYZ = val;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasXpos(G4double val)
{
  fXposGas  = val;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetWorldSizeX(G4double val)
{
  fWorldSizeX = val;
  fDefaultWorld = false;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  fWorldSizeYZ = val;
  fDefaultWorld = false;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetPairEnergy(G4double val)
{
  //Debugging
  // G4cout << fGasMaterial << G4endl;
  if(val > 0.0) {
    fGasMaterial->GetIonisation()->SetMeanEnergyPerIonPair(val);
  }
  //G4cout << fGasMaterial << G4endl;
}

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Functions to set anode geometry
void DetectorConstruction::SetAnodePosition(G4double position)
{
  fAnodePosition = position;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAnodeLength(G4double length)
{
  fAnodeX = length;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSegmentLength(G4double SegLength)
{
  fSegmentX = SegLength;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
