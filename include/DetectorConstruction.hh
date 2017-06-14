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
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4Region.hh" //Required for PAI model
#include "G4ProductionCuts.hh" //Required for PAI model
#include "TargetSD.hh"

class G4Box;
class G4Material;
class G4VPhysicalVolume;
class ostringstream;
class G4MaterialCutsCouple;
class F02ElectricFieldSetup;
class DetectorMessenger;
class G4NistManager;
//class SegParameterisation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();
public:

     G4Material* GetMaterial() {return fGasMaterial;};

  public:

     void SetWindowMaterial (G4String);
     void SetWindowThickness(G4double);

  // void SetGasMaterial (G4Material*);
     void SetGasMaterial(G4String);
     void SetGasThickness(G4double);
     void SetGasSizeYZ   (G4double);

     void SetGasXpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeX   (G4double);
     void SetWorldSizeYZ  (G4double);

     void SetAnodePosition(G4double);
     void SetAnodeLength  (G4double);
     void SetSegmentLength(G4double);
  // void SetMagField(G4double);

     virtual G4VPhysicalVolume* Construct();
  //virtual void ConstructField();
     void SetPairEnergy(G4double);
     void UpdateGeometry();

     void SetGasPressure(G4double);

     void PrintCalorParameters();

     G4Material* GetWindowMaterial()	{return fWindowMaterial;};
     G4double    GetWindowThickness()   {return fWindowThickness;};
     G4double    GetWindowSizeYZ()	{return fWindowSizeYZ;};

     G4Material* GetGasMaterial()  {return fGasMaterial;};
     G4double    GetGasThickness() {return fGasThickness;};
     G4double    GetGasSizeYZ()    {return fGasSizeYZ;};

     G4double    GetGasXpos()      {return fXposGas;};
     G4double    GetxstartGas()         {return fXstartGas;};
     G4double    GetxendGas()           {return fXendGas;};

     G4double    GetNumberOfSegments()     {return nSegments;};
     G4double    GetSegmentThickness()     {return fSegmentX;};

     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeX()        {return fWorldSizeX;};
     G4double    GetWorldSizeYZ()       {return fWorldSizeYZ;};

     G4double    GetDSSSDDetectorSizeYZ()            {return fDSSSDDetectorSizeYZ;};

     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};
     const G4VPhysicalVolume* GetGas()   {return fPhysiGas;};
     const G4VPhysicalVolume* GetSegment() {return fPhysiSegment;};
  // const G4VPVParameterisation* GetSegment() { return fPhysiParam;};
     const G4VPhysicalVolume* GetWindow()     {return fPhysiWindow;};
     const G4VPhysicalVolume* GetDSSSD()      {return fPhysiDSSSDActive;};
     const G4MaterialCutsCouple* GetGasMaterialCut()  const
                             {return fLogicGas->GetMaterialCutsCouple();};
  private:
    G4NistManager* man;

  // G4double IsobutaneDensity, IsobutanePressure;
     //DSSSD Components
     // active silicon strip
     G4Box*             fSolidDSSSDActive;
     G4LogicalVolume*   fLogicDSSSDActive;
     G4VPhysicalVolume* fPhysiDSSSDActive;
     G4Material*	fDSSSDActiveMaterial;
     G4double		fDSSSDActiveSizeZ;	// active strip width
     G4double		fDSSSDActiveThickness;	// active, gap, holder same thickness

     // dead layer
     G4Box*             fSolidDSSSDDeadlayer;
     G4LogicalVolume*   fLogicDSSSDDeadlayer;
     G4VPhysicalVolume* fPhysiDSSSDDeadlayer;
     G4Material*	fDSSSDDeadlayerMaterial;
     G4double		fDSSSDDeadlayerThickness;

     // mother holder for front and backs strips and deadlayer
     G4Box*             fSolidDSSSDDetector;
     G4LogicalVolume*   fLogicDSSSDDetector;
     G4VPhysicalVolume* fPhysiDSSSDDetector;
     G4Material*	fDSSSDDetectorMaterial;
     G4double		fDSSSDDetectorThickness; // back layer + front layer + dead layer
     G4double		fDSSSDDetectorSizeYZ;
     G4double		fXposDSSSDDetector;

/*   BLOCKED UNTIL POSITION NEEDED
     // gap between strips
     G4Box*             fSolidDSSSDGap;
     G4LogicalVolume*   fLogicDSSSDGap;
     G4VPhysicalVolume* fPhysiDSSSDGap;
     G4Material*	fDSSSDGapMaterial;
     G4double		fDSSSDGapSizeZ;		// gap width
     // mother holder for Si + gap and then gets replicated for front and back strips
     G4Box*             fSolidDSSSDStrip;
     G4LogicalVolume*   fLogicDSSSDStrip;
     G4VPhysicalVolume* fPhysiDSSSDStrip;
     G4Material*	fDSSSDStripMaterial;
     G4double		fDSSSDStripThickness;	// active, gap, holder same thickness
     G4double		fDSSSDStripSizeY;	// strip length
     G4double		fDSSSDStripSizeZ;	// strip width
      // mother holder for front strips
     G4Box*             fSolidDSSSDFront;
     G4LogicalVolume*   fLogicDSSSDFront;
     G4VPhysicalVolume* fPhysiDSSSDFront;
     G4Material*	fDSSSDFrontMaterial;
     // mother holder for back strips
     G4Box*             fSolidDSSSDBack;
     G4LogicalVolume*   fLogicDSSSDBack;
     G4VPhysicalVolume* fPhysiDSSSDBack;
     G4Material*	fDSSSDBackMaterial;
*/
     //Chamber Window
     G4Box*             fSolidWindow;
     G4LogicalVolume*   fLogicWindow;
     G4VPhysicalVolume* fPhysiWindow;
     G4Material*	fWindowMaterial;
     G4double       fWindowThickness;
     G4double		fWindowSizeYZ;
     G4double		fXposWindow;


     // Ionization Gas Gas
     G4Box*             fSolidGas;
     G4LogicalVolume*   fLogicGas;
     G4VPhysicalVolume* fPhysiGas;
     G4Material*        fGasMaterial;
     G4double           fGasThickness;
     G4double           fGasSizeYZ;
     G4double           fXposGas;
     G4double           fXstartGas, fXendGas;


  //Anode Holder Volume
  G4double fAnodeX;
  G4double fAnodePosition;
  G4Box* fSolidAnodeBox;
  G4LogicalVolume* fLogicAnodeBox;
  G4VPhysicalVolume* fPhysiAnodeBox;

  //Anode Segmentation
  G4double fSegmentX;
  G4int nSegments;
  G4double fAnodeDeadLayerThickness;
  G4Box* fSolidSegment;
  G4LogicalVolume* fLogicSegment;
  G4VPhysicalVolume* fPhysiSegment;


     // World
     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;
     G4Material*        fWorldMaterial;
     G4Region* fRegGasDet; //Required for PAI model
     G4ProductionCuts* fGasDetectorCuts; //Required for PAI model
     TargetSD*          fTargetSD;
     G4double           fWorldSizeX;
     G4double           fWorldSizeYZ;
     G4bool             fDefaultWorld;

  //G4UniformMagField* fMagField;
  G4Cache<F02ElectricFieldSetup*> fEmFieldSetup;
     DetectorMessenger* fDetectorMessenger;

  G4bool fCheckOverlaps; //Option to activate checking of overlaps
  private:

     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
