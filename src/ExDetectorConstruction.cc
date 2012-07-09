//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExDetectorConstruction.cc,v 1.8 2003/10/06 08:59:11 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include "ExDetectorConstruction.hh"
#include "ExTrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UnionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbVolumeMgr.hh"
#include "ExGlobalParameters.hh"
#include "globals.hh"

ExDetectorConstruction::ExDetectorConstruction()
: experimentalHall_log(0), CZTDet_log(0), CZTDead_log(0),
  experimentalHall_phys(0), CZTDet_phys(0), CZTDead_phys(0),
  collTungstenLogic(0), collTungstenPhys(0),
  collLeadLogic(0), collLeadPhys(0),
  collAluminumBigLogic(0), collAluminumBigPhys(0),
  collAluminumSmallLogic(0), collAluminumSmallPhys(0),
  collAluminumCoverLogic(0), collAluminumCoverPhys(0),
  holderUnionVolumeLogic(0), holderUnionVolumePhys(0), 
  sourceHolderExtensionLogic(0), sourceHolderExtensionPhys(0),
  firstSubtractLogic(0),firstSubtractPhys(0),
  sourceHolderCoverLogic(0), sourceHolderCoverPhys(0),
  sourceHolderCo57LeadLogic(0), sourceHolderCo57LeadPhys(0),
  AlBoxCoverLogic(0), AlBoxCoverPhys(0),
  AlCs137HolderLogic(0), AlCs137HolderPhys(0),
  sourceHolderExtensionCs137Logic(0), sourceHolderExtensionCs137Phys(0),
  LeadInsideCs137Logic(0), LeadInsideCs137Phys(0),
  sourceHolderCoverCs137Logic(0), sourceHolderCoverCs137Phys(0),
  collimatorTestLogic(0), collimatorTestPhys(0)
{
  messenger = new G4tgrMessenger;
}
  
ExDetectorConstruction::~ExDetectorConstruction()
{
  delete messenger;
}
  
G4VPhysicalVolume* ExDetectorConstruction::Construct()
{
  
  // ------------------------------------------------------ materials 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  
  // n=number of nucleons in an isotope;
  G4int ncomponents, natoms;
  G4double fractionmass;
  
  // Collimator and shielding materials.
  // Pb and W and Fe and Al.
  G4Material* Pb = new G4Material("Lead", z= 82., a= 207.19*g/mole, density= 11.35*g/cm3);
  G4Material* W = new G4Material("Tungsten", z= 74., a= 183.84*g/mole, density= 19.25*g/cm3);
  G4Material* Fe = new G4Material("Iron", z= 26., a= 55.845*g/mole, density= 7.874*g/cm3);
  G4Material* Al =new G4Material("Aluminum", z=13., a= 27.*g/mole, density= 2.7*g/cm3 );

  // Air and H2O
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  
  G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  
  // Vacuum
  G4Material* Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density, kStateGas, 3.e-18*pascal, 2.73*kelvin);
  
  // ------------------------------------------------------ //
  // Detector materials.
  // CdZnTe
  G4Element* Cd = new G4Element("Cadmium", symbol="Cd", z=48., a=112.411*g/mole);
  G4Element* Zn = new G4Element("Zinc", symbol="Zn", z=30., a=65.39*g/mole);
  G4Element* Te = new G4Element("Tellurium", symbol="Te", z=52., a=127.6*g/mole);
  
  G4Material* CZT = new G4Material("CZT", density=5.78*g/cm3, ncomponents=3);
  CZT->AddElement(Cd, natoms=90);
  CZT->AddElement(Zn, natoms=10);
  CZT->AddElement(Te, natoms=100);
  
  // HgI2
  G4Element* Hg = new G4Element("Mercury", symbol="Hg", z=80., a=200.59*g/mole);
  G4Element* I = new G4Element("Iodine", symbol="I", z=53., a=126.90447*g/mole);
  
  G4Material* HgI2 = new G4Material("HgI2", density=6.4*g/cm3, ncomponents=2);
  HgI2->AddElement(Hg, natoms=1);
  HgI2->AddElement(I, natoms=2);
  
  // TlBr. Added by B.D. 14/09/10.
  G4Element* Tl = new G4Element("Thallium", symbol="Tl", z=81., a=204.38*g/mole);
  G4Element* Br = new G4Element("Bromine", symbol="Br", z=35., a=79.9*g/mole);
  
  G4Material* TlBr = new G4Material("TlBr", density=7.58*g/cm3, ncomponents=2);
  TlBr->AddElement(Tl, natoms=1);
  TlBr->AddElement(Br, natoms=1);
  // Detector materials ends here!
  // ------------------------------------------------------ //
  

  // ------------------------------------------------------ volumes


  // ------------------------------ experimental hall (world volume)  
  // I never use this, so I don't know what exactly it does.
  // I kept it as it is just in case. Note added by B.D. 21/09/10.
  if (CGlobalParameters::enableTextGeom)
  {
    // if use text geometry to define the environment
    G4String filename = "g4geom.txt";
    G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
    volmgr->AddTextFile(filename);
    experimentalHall_phys = volmgr->ReadAndConstructDetector();
    experimentalHall_log = experimentalHall_phys->GetLogicalVolume();
  }
  else
  {
    // World volume is created using the following setup.
    // for simple use
    G4double expHall_x = 3.0*m;
    G4double expHall_y = 3.0*m;
    G4double expHall_z = 3.0*m;
    G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
    experimentalHall_log = new G4LogicalVolume(experimentalHall_box, Vacuum,"expHall_log",0,0,0);
    experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(), experimentalHall_log,"expHall",0,false,0);
    experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
   }
  // World volume ends here!
  // ------------------------------------------------------ //

  
  // -------------- Detector Definitions! ---------------------//
  // Detector size.
  G4double det_x = CGlobalParameters::detsizeX;
  G4double det_y = CGlobalParameters::detsizeY;
  G4double det_z = CGlobalParameters::detsizeZ;
  
  // Dead layer definition in z. Comment added by B.D. 14/09/10.
  // Might be useful in the future, so I kept it.
  G4double det_dead = CGlobalParameters::detdeadZ;
  det_z -= det_dead;
  
  // Below is to create array of detectors.
  // Might be useful in the future, so I kept it.
  G4double detintv_x = CGlobalParameters::detintvX;
  G4double detintv_y = CGlobalParameters::detintvY;
  G4double detintv_z = CGlobalParameters::detintvZ;
  
  G4int detmatrix_x = (G4int)CGlobalParameters::detmatrixX;
  G4int detmatrix_y = (G4int)CGlobalParameters::detmatrixY;
  G4int detmatrix_z = (G4int)CGlobalParameters::detmatrixZ;
  
  // Type of the detector. CZT, HgI2 or TlBr.
  G4int dettype = (G4int)CGlobalParameters::detType;
  
  // ------------------------------------------------------ //
  // Detector geometry and layout in the world volume.
  // Start!AlBoxCoverLogic
  G4Box* CZTDet_box = new G4Box("CZTDet_box", det_x/2, det_y/2, det_z/2);
  
  G4Box* CZTDead_box = 0;
  if (det_dead>0)
  {
    CZTDead_box = new G4Box("CZTDead_box", det_x/2, det_y/2, det_dead/2);
  }
  
  switch (dettype)
  { 
    case 0:
      CZTDet_log = new G4LogicalVolume(CZTDet_box, CZT,"CZTDet_log",0,0,0);
      if (det_dead>0)
      {
        CZTDead_log = new G4LogicalVolume(CZTDead_box, CZT,"CZTDead_log",0,0,0);
      }
      break;
    case 1:
      CZTDet_log = new G4LogicalVolume(CZTDet_box, HgI2,"CZTDet_log",0,0,0);
      if (det_dead>0)
      {
        CZTDead_log = new G4LogicalVolume(CZTDead_box, HgI2,"CZTDead_log",0,0,0);
      }
      break;
    case 2:
      CZTDet_log = new G4LogicalVolume(CZTDet_box, TlBr,"CZTDet_log",0,0,0);
      if (det_dead>0)
      {
        CZTDead_log = new G4LogicalVolume(CZTDead_box, TlBr,"CZTDead_log",0,0,0);
      }
      break;
  }

  // G4AssemblyVolume for the array
  // Define detector array as assembly volumes
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  
  // Rotation and translation of a plate inside the assembly
  G4RotationMatrix * Ra = new G4RotationMatrix();
  G4ThreeVector Ta;
  
  // Rotation of the assembly inside the world
  G4RotationMatrix * Rm = new G4RotationMatrix();
  G4ThreeVector Tm;
  
  // Fill the assembly by the Detector
  Ta = G4ThreeVector(0,0,0);
  assemblyDetector->AddPlacedVolume( CZTDet_log, Ta, Ra);
  if (det_dead != 0)
  {
    Ta = G4ThreeVector(0, 0, -(det_dead + det_z)/2);
    assemblyDetector->AddPlacedVolume(CZTDead_log, Ta, Ra);
  }
  
  for (int i=0; i<detmatrix_x; i++)
  {
    for (int j=0; j<detmatrix_y; j++)
    {
      for (int k=0; k<detmatrix_z; k++)
      {
        Tm.setX(detintv_x*(i-(detmatrix_x-1)/2.0));
        Tm.setY(detintv_y*(j-(detmatrix_y-1)/2.0));
        Tm.setZ((detintv_z+(det_z+det_dead))*(k-(detmatrix_z-1)/2.0));

        assemblyDetector->MakeImprint(experimentalHall_log, Tm, Rm); 
      }
    }
  }
  // End!
  // ------------------------------------------------------ //


  // ******************************************************************** //
  // Collimator geometry added. O.A. and B.D. 14/09/10.
  // ******************************************************************** //

  // Needed to shift geometries. O.A. and B.D. 15/09/10.
  G4double shift;
  G4double shiftCollimator; // Needed to shift holder.
  // Reference point to shift collimator and source holder with respect to the detector.
  G4double ReferenceFrame;
  
  // Define Al Box cover distance.
  G4double AlBoxCoverDistance = CGlobalParameters::alboxcoverdistance;
  
  // Collimator and/or source holder's distance to the detector
  // Unit conversion added. There seems to be some problem when we add a ReferenceFrame.
  // This solved the problem.
  G4double distanceToTheDetector = CGlobalParameters::distancetothedetector+AlBoxCoverDistance;
  
  // ******************************************************************** //
  // Al box cover geometry added. O.A. and B.D. 20/09/10.
  // It can be turn on and off from settings.ini file.
  // ******************************************************************** //
  G4int enableAlBox = (G4int)CGlobalParameters::enableAlBoxCover;
  
  if (enableAlBox)
  {
    // Al box cover dimension.
    G4double AlBox_x = 9.0*cm/2.0;
    G4double AlBox_y = 11.5*cm/2.0;
    G4double AlBox_z = 3.0*mm/2.0;
    
    // Cut out piece.
    G4double CutAlBox_x = 3.5*cm/2.0;
    G4double CutAlBox_y = 4.0*cm/2.0;
    G4double CutAlBox_z = 3.01*mm/2.0;
    
    // Al box cover and cut out piece definition.
    G4Box* AlBoxCover = new G4Box("AlBoxCover",AlBox_x,AlBox_y,AlBox_z);
    G4Box* CutAlBoxCover = new G4Box("CutAlBoxCover",CutAlBox_x,CutAlBox_y,CutAlBox_z);

    // Actually there is no rotation. It is needed for solid subtraction.
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateZ(0*deg);
    zRot->rotateX(0*deg);
    zRot->rotateY(0*deg);
    
    // Cut out piece 0.6cm inside. 
    G4ThreeVector zTrans(0,CutAlBox_y-AlBox_y+0.6*cm,0);

    G4SubtractionSolid*  AlBoxCoverVolume = new G4SubtractionSolid("AlBoxCoverVolume", AlBoxCover, CutAlBoxCover, zRot, zTrans);
 
    // Position of Al box cover.
    G4ThreeVector AlBoxCoverPos = G4ThreeVector(0,3.1*cm,AlBoxCoverDistance);
    
    // Definition of logical and physical volume.
    AlBoxCoverLogic = new G4LogicalVolume(AlBoxCoverVolume, Al,"AlBoxCoverLogic",0,0,0);
    AlBoxCoverPhys = new G4PVPlacement(0, AlBoxCoverPos,AlBoxCoverLogic,"AlBoxCoverPhys",experimentalHall_log,false,0);
    
    // Is volume visible?
    G4VisAttributes* visibleAlBoxCover = new G4VisAttributes(G4Colour(1,0,0)); 
    AlBoxCoverLogic->SetVisAttributes(visibleAlBoxCover);
//   AlBoxCoverLogic->SetVisAttributes(G4VisAttributes::Invisible);
  }  
  // ******************************************************************** //
  // Al box cover ends here!
  // ******************************************************************** //
  
  
  // ******************************************************************** //
  // Collimator geometry added. O.A. and B.D. 15/09/10.
  // It can be turn on and off from settings.ini file.
  // ******************************************************************** //
  

   if (CGlobalParameters::enableCollimator == 0)
  {
    G4cout << " -----------------------------------------------------------" << G4endl;
    G4cout << "<<<<---- Collimator is NOT present! ---->>>>" << G4endl;
  }    

  if (CGlobalParameters::enableCollimator == 1 )
  {
    G4cout << " -----------------------------------------------------------" << G4endl;
    G4cout << "<<<<---- Collimator is present! ---->>>>" << G4endl;
    
    // Aluminum Cover towards the detector, i.e. tip of the collimator.
    // This piece taken as a reference.
    // All collimator and source holder parts will be shifted according to this piece.
    
    // Define dimensions.
    G4double innerRadiusOfTheTubeAlC = 0.7*mm;
    G4double outerRadiusOfTheTubeAlC = 22.5*mm; 
    G4double heightOfTheTubeAlC = 0.8*cm/2.0;
    G4double startAngleOfTheTube = 0*deg;
    G4double spanningAngleOfTheTube = 360*deg;
    
    // Position of the cover.
    ReferenceFrame = heightOfTheTubeAlC+det_z/2.0;
    G4double coverAlX = 0.0*cm;
    G4double coverAlY = 0.0*cm;
    G4double coverAlZ = distanceToTheDetector+ReferenceFrame;
    
    // Define the volume. 
    G4Tubs* collAluminumCover = new G4Tubs("Coll_AluminumCover", innerRadiusOfTheTubeAlC,
        outerRadiusOfTheTubeAlC, heightOfTheTubeAlC, startAngleOfTheTube, spanningAngleOfTheTube);
    
    // Position vector.
    G4ThreeVector collAluminumCoverPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ);
   
    // Definition of logical and physical volume.    
    collAluminumCoverLogic = new G4LogicalVolume(collAluminumCover, Al, "collAluminumCover",0,0,0);
    collAluminumCoverPhys = new G4PVPlacement(0, collAluminumCoverPos, collAluminumCoverLogic,
        "collAluminumCoverPhys", experimentalHall_log, false, 0);
 
    // Is volume visible?    
    G4VisAttributes* visibleAluminumCover = new G4VisAttributes(G4Colour(0.0,1.0,1.0)); 
    collAluminumCoverLogic->SetVisAttributes(visibleAluminumCover);
    //collAluminumCoverLogic->SetVisAttributes(G4VisAttributes::Invisible);
    // ------------------------------------------------------ //
    // Tungsten part of the collimator.
    // Define dimensions.
    G4double innerRadiusOfTheTube = 0.35*mm;
    G4double outerRadiusOfTheTube = 1.15*mm;
    G4double heightOfTheTube = 0.5*cm/2.0;

    // Define the volume. 
    G4Tubs* collTungsten = new G4Tubs("Coll_Tungsten", innerRadiusOfTheTube,
        outerRadiusOfTheTube, heightOfTheTube, startAngleOfTheTube, spanningAngleOfTheTube);

    // Position of tungsten.    
    G4ThreeVector collTungstenPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ+heightOfTheTube+heightOfTheTubeAlC);    

    // Definition of logical and physical volume.    
    collTungstenLogic = new G4LogicalVolume(collTungsten, W, "collTungsten",0,0,0);
    collTungstenPhys = new G4PVPlacement(0, collTungstenPos, collTungstenLogic,
        "collTungstenPhys", experimentalHall_log, false, 0);
   
    // Is volume visible?        
    G4VisAttributes* visibleTungsten = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); 
    collTungstenLogic->SetVisAttributes(visibleTungsten);
    //collTungstenLogic->SetVisAttributes(G4VisAttributes::Invisible);
    // ------------------------------------------------------ //
    
    // ------------------------------------------------------ //
    // Lead part of the collimator.
    // Define dimensions.
    G4double innerRadiusOfTheTubeL = 1.15*mm;
    G4double outerRadiusOfTheTubeL = 10.0*mm;
    
    // Define the volume.
    G4Tubs* collLead = new G4Tubs("Coll_Lead", innerRadiusOfTheTubeL,
        outerRadiusOfTheTubeL, heightOfTheTube, startAngleOfTheTube, spanningAngleOfTheTube);
    
    // Position of lead.
    G4ThreeVector collLeadPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ+heightOfTheTube+heightOfTheTubeAlC);
    
    // Definition of logical and physical volume.        
    collLeadLogic = new G4LogicalVolume(collLead, Pb, "collLead",0,0,0);
    collLeadPhys = new G4PVPlacement(0, collLeadPos, collLeadLogic,
        "collLeadPhys", experimentalHall_log, false, 0);

    // Is volume visible?            
    G4VisAttributes* visibleLead = new G4VisAttributes(G4Colour(144,144,144)); 
    collLeadLogic->SetVisAttributes(visibleLead);
    //collLeadLogic->SetVisAttributes(G4VisAttributes::Invisible);
    // ------------------------------------------------------ //
    
    // ------------------------------------------------------ //
    // Aluminum Cover.
    // Covers lead and tungsten parts.
    // Define dimensions.
    G4double innerRadiusOfTheTubeAlB = 10.0*mm;
    G4double outerRadiusOfTheTubeAlB = 22.5*mm;
    G4double heightOfTheTubeAlB = 0.5*cm/2.0;
    
    // Define the volume.
    G4Tubs* collAluminumBig = new G4Tubs("Coll_AluminumBig", innerRadiusOfTheTubeAlB,
        outerRadiusOfTheTubeAlB, heightOfTheTubeAlB, startAngleOfTheTube, spanningAngleOfTheTube);
    
    // Position of the cover.    
    G4ThreeVector collAluminumBigPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ+heightOfTheTubeAlB+heightOfTheTubeAlC);
    
    // Definition of logical and physical volume.        
    collAluminumBigLogic = new G4LogicalVolume(collAluminumBig, Al, "collAluminumBig",0,0,0);
    collAluminumBigPhys = new G4PVPlacement(0, collAluminumBigPos, collAluminumBigLogic,
        "collAluminumBigPhys", experimentalHall_log, false, 0);
    
    // Is volume visible?            
    G4VisAttributes* visibleAluminumBig = new G4VisAttributes(G4Colour(248,248,248)); 
    collAluminumBigLogic->SetVisAttributes(visibleAluminumBig);
    //collAluminumBigLogic->SetVisAttributes(G4VisAttributes::Invisible);
    // ------------------------------------------------------ //
    
    // ------------------------------------------------------ //
    // Aluminum cap.
    G4double innerRadiusOfTheTubeAlS = 1.0*mm;
    G4double outerRadiusOfTheTubeAlS = 22.5*mm;
    G4double heightOfTheTubeAlS = 1.2*cm/2.0;
    shift = ReferenceFrame+heightOfTheTubeAlB+heightOfTheTubeAlC+2*mm;
      
    // Define the volumes.
    G4Tubs* subtractVolume = new G4Tubs("subtractVolume",0.75*cm,outerRadiusOfTheTubeAlS+0.1*mm,1.0*cm/2.0,
        startAngleOfTheTube, spanningAngleOfTheTube);
    G4Tubs* collAluminumSmallBefore = new G4Tubs("Coll_AluminumSmallBefore", innerRadiusOfTheTubeAlS,
        outerRadiusOfTheTubeAlS, heightOfTheTubeAlS, startAngleOfTheTube, spanningAngleOfTheTube);

    // Actually there is no rotation. It is needed for solid subtraction.    
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateZ(0*deg);
    zRot->rotateX(0*deg);
    zRot->rotateY(0*deg);
    // Needed for volume subtraction. 
    G4ThreeVector zTrans(0,0,+0.2*cm/2.0);
    
    // This is the remaining volume to be used.
    G4SubtractionSolid* collAluminumSmall = new G4SubtractionSolid("Coll_AluminumSmall", collAluminumSmallBefore, subtractVolume, zRot, zTrans);
    
    // Position of the aluminum cap.
    G4ThreeVector collAluminumSmallPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ+shift);

    // Definition of logical and physical volume.        
    collAluminumSmallLogic = new G4LogicalVolume(collAluminumSmall, Al, "collAluminumSmall",0,0,0);   
    collAluminumSmallPhys = new G4PVPlacement(0, collAluminumSmallPos, collAluminumSmallLogic,
        "collAluminumSmallPhys", experimentalHall_log, false, 0);
    
    // Is volume visible?            
    G4VisAttributes* visibleAluminumSmall = new G4VisAttributes(G4Colour(248,248,248)); 
    collAluminumSmallLogic->SetVisAttributes(visibleAluminumSmall);
   // collAluminumSmallLogic->SetVisAttributes(G4VisAttributes::Invisible);
    // ------------------------------------------------------ //
  }  

  if (CGlobalParameters::enableCollimator == 2)
  {
    G4cout << " -----------------------------------------------------------" << G4endl;
    G4cout << "<<<<---- Collimator is a test collimator! ---->>>>" << G4endl;

     // Lead for collimator test.
     // Define dimensions.
    G4double innerRadiusOfTheTest = CGlobalParameters::collimatorTestHoleDiameter/2.0;
    G4double outerRadiusOfTheTest = CGlobalParameters::collimatorTestDiameter/2.0;
    G4double heightOfTheTest = CGlobalParameters::collimatorTestThickness/2.0;
    G4double startAngleOfTheTest = 0*deg;
    G4double spanningAngleOfTheTest = 360*deg;
    
     // Position of the cover.
    ReferenceFrame = heightOfTheTest+det_z/2.0;
    G4double coverAlX = 0.0*cm;
    G4double coverAlY = 0.0*cm;
    G4double coverAlZ = distanceToTheDetector+ReferenceFrame;
    
    // Define volumes.
    G4Tubs* collimatorTest = new G4Tubs("collimatorTest",
        innerRadiusOfTheTest,
        outerRadiusOfTheTest,
        heightOfTheTest,
        startAngleOfTheTest,
        spanningAngleOfTheTest);
 
    // Position of lead.
    G4ThreeVector collimatorTestPos = G4ThreeVector(coverAlX, coverAlY, coverAlZ);
    
    // Definition of logical and physical volume.        
    collimatorTestLogic = new G4LogicalVolume(collimatorTest, Pb, "collimatorTest",0,0,0);
    collimatorTestPhys = new G4PVPlacement(0, collimatorTestPos, collimatorTestLogic,
        "collimatorTestPhys", experimentalHall_log, false, 0);
    
    // Is volume visible?            
    G4VisAttributes* visibleTest = new G4VisAttributes(G4Colour(144,144,144)); 
    collimatorTestLogic->SetVisAttributes(visibleTest);
 
  }  
  
  // ******************************************************************** //
  // Collimator definition ends here!
  // ******************************************************************** //
  

  // ******************************************************************** //
  // Source Holder geometry added. O.A. and B.D. 15/09/10.
  // ******************************************************************** //

  G4int sourceHolderType = 2;//(G4int)CGlobalParameters::enableSourceHolder;
  
  switch (sourceHolderType)
  {
    case 0: // There is no source holder.
       G4cout << "<<<<---- No source holder. ---->>>>" << G4endl;
       G4cout << " -----------------------------------------------------------" << G4endl;
       break;

    case 1: // Source holder geometry is for Am, Co or Cd.
       {
         G4cout << "Source holder for Am, Co or Cd." << G4endl;
         G4cout << " -----------------------------------------------------------" << G4endl;
         
         // ------------------------------------------------------ //
         // Define dimensions.
         G4double innerRadiusOfTheHolderFirst = 0.2*cm;
         G4double outerRadiusOfTheHolder= 22.5*mm;
         G4double heightOfTheHolderFirst = 4.5*cm/2.0;
         G4double startAngleOfTheHolder = 0*deg;
         G4double spanningAngleOfTheHolder = 360*deg;
         
         // This is the part that will be subtracted.
         G4double innerRadiusOfTheHolderSecond = 0.0*mm;
         G4double outerRadiusOfTheHolderSecond = 1.75*cm;
         G4double heightOfTheHolderSecond = 3.0005*cm/2.0;
         
         ReferenceFrame = det_z/2.0+heightOfTheHolderFirst;
         
         G4double holderAlX = 0.0*cm;
         G4double holderAlY = 0.0*cm; 
         G4double holderAlZ = distanceToTheDetector+ReferenceFrame;
 
         // ShiftHolder. Needed below!.
         // 1.5*cm is the thickness of the collimator for Am, Co and Cd.
         // Shift depends on the existence of collimator.
         if ((CGlobalParameters::enableCollimator == 0) || (CGlobalParameters::enableCollimator == 2)) 
         {
           shiftCollimator = 1.0*cm;
         } else if (CGlobalParameters::enableCollimator == 1)
         {
           shiftCollimator = 2.5*cm;
         } 
         
         // Define volumes.
         G4Tubs* firstUnionVolume = new G4Tubs("firstUnionVolume",
             innerRadiusOfTheHolderFirst,
             outerRadiusOfTheHolder,
             heightOfTheHolderFirst,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* secondUnionVolume = new G4Tubs("secondUnionVolume",
             innerRadiusOfTheHolderSecond,
             outerRadiusOfTheHolderSecond,
             heightOfTheHolderSecond,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         // Actually there is no rotation. It is needed for solid subtraction.
         G4RotationMatrix* zRot = new G4RotationMatrix;
         zRot->rotateZ(0*deg);
         zRot->rotateX(0*deg);
         zRot->rotateY(0*deg);
         
         // This is the shift needed to subtract correct part of the volume.
         G4ThreeVector zTrans(0,0,heightOfTheHolderSecond/2.0);
    
         // This is the remaining volume to be used.
         G4SubtractionSolid*  holderUnionVolume = new G4SubtractionSolid("holderUnionVolume", firstUnionVolume, secondUnionVolume, zRot, zTrans);
         
         // Position of the holder.
         G4ThreeVector holderUnionVolumePos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator);
         
         // Define logical and physical volume.
         holderUnionVolumeLogic = new G4LogicalVolume(holderUnionVolume, Al, "holderUnionVolume",0,0,0);
         holderUnionVolumePhys = new G4PVPlacement(0, holderUnionVolumePos, holderUnionVolumeLogic,
             "holderUnionVolumePhys", experimentalHall_log, false, 0);
         
         // Is volume visible?            
         G4VisAttributes* visibleholderUnionVolume = new G4VisAttributes(G4Colour(248,248,248)); 
         holderUnionVolumeLogic->SetVisAttributes(visibleholderUnionVolume);
         // holderUnionVolumeLogic->SetVisAttributes(G4VisAttributes::Invisible);//
         // ------------------------------------------------------ //
         
         // ------------------------------------------------------ //
         // Lead inside the source holder.
         // Define dimensions.
         G4double innerRadiusOfTheLeadInside1 = 0.0*mm;
         G4double outerRadiusOfTheLeadInside = 1.75*cm;
         G4double heightOfTheLeadInside1 = 1.984*cm/2.0;

         G4double outerRadiusOfTheLeadInside2 = 0.2*cm;
         G4double heightOfTheLeadInside2 = 0.992*cm/2.0;
         
         G4double outerRadiusOfTheLeadInside3 = 3.96*mm; // (For the old source 1.25 cm, for new source 3.96*mm)
         G4double heightOfTheLeadInside3 = 0.496*cm/2.0; 

         // Define volumes.
         G4Tubs* firstLeadInside = new G4Tubs("firstLeadInside",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside,
             heightOfTheLeadInside1,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);

         G4Tubs* secondLeadInside = new G4Tubs("secondLeadInside",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside2,
             heightOfTheLeadInside2,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);

         G4Tubs* thirdLeadInside = new G4Tubs("thirdLeadInside",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside3,
             heightOfTheLeadInside3,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
     
         // Subtract second lead from the first.
         G4ThreeVector zTransLead1(0,0,-heightOfTheLeadInside2);
         G4SubtractionSolid*  firstSubtractVolume = new G4SubtractionSolid("firstSubtractVolume", firstLeadInside, secondLeadInside, zRot, zTransLead1);

         // Subtract third lead from the remaining.
         G4ThreeVector zTransLead2(0,0,+heightOfTheLeadInside3/2.0);
         G4SubtractionSolid*  secondSubtractVolume = new G4SubtractionSolid("secondSubtractVolume", firstSubtractVolume, thirdLeadInside, zRot, zTransLead2);
         
         // Position of the holder.         
         G4ThreeVector firstSubtractPos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator+heightOfTheLeadInside3);
         
         // Definition of logical and physical volume.
         firstSubtractLogic = new G4LogicalVolume(secondSubtractVolume, Pb, "firstSubtractVolume",0,0,0);
         firstSubtractPhys = new G4PVPlacement(0, firstSubtractPos, firstSubtractLogic,
             "firstSubtractPhys", experimentalHall_log, false, 0);
                
         // Is volume visible?          
         G4VisAttributes* visiblefirstSubtract = new G4VisAttributes(G4Colour(144,144,144)); 
         firstSubtractLogic->SetVisAttributes(visiblefirstSubtract);
	//firstSubtractLogic->SetVisAttributes(G4VisAttributes::Invisible);
         // ------------------------------------------------------ //

         
         // ------------------------------------------------------ //
         // Cap of the source holder.
         // Define dimensions.
         G4double innerRadiusOfTheSourceHolderCover1 = 1.75*cm;
         G4double outerRadiusOfTheSourceHolderCover1 = 2.344*cm;
         G4double heightOfTheSourceHolderCover1 = 1.19*cm/2.0;

         G4double innerRadiusOfTheSourceHolderCover2 = 0.0*mm;
         G4double outerRadiusOfTheSourceHolderCover2 = 2.304*cm;
         G4double heightOfTheSourceHolderCover2 = 1.84*cm/2.0;
         
         // Define volumes.
         G4Tubs* firstSourceHolderSub = new G4Tubs("firstSourceHolderSub",
             innerRadiusOfTheSourceHolderCover1,
             outerRadiusOfTheSourceHolderCover1,
             heightOfTheSourceHolderCover1,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* firstSourceHolderCover = new G4Tubs("firstSourceHolderCover",
             innerRadiusOfTheSourceHolderCover2,
             outerRadiusOfTheSourceHolderCover2,
             heightOfTheSourceHolderCover2,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);

         // This is the shift needed to subtract correct part of the volume.
         G4ThreeVector zTransCover(0,0,+heightOfTheSourceHolderCover1-heightOfTheSourceHolderCover2-0.01*mm);
    
         // This is the volume to be used.
         G4SubtractionSolid* sourceHolderCover = new G4SubtractionSolid("sourceHolderCover", firstSourceHolderCover, firstSourceHolderSub, zRot, zTransCover);
         
         // Define logical volume of the cover.
         sourceHolderCoverLogic = new G4LogicalVolume(sourceHolderCover, Al, "sourceHolderCover",0,0,0);
         
         // Source holder geometry is for Co 57.
         // Add one more lead piece if the source holder is for Co.
         if (CGlobalParameters::sourceHolderIsCo57)
         {
           
           // Lead piece for Co 57 closer to the source holder cap.
           // Define dimensions.
           G4double innerRadiusOfTheLeadCo57 = 0.0*mm;
           G4double outerRadiusOfTheLeadCo57 = 1.75*cm;
           G4double heightOfTheLeadCo57 = 0.496*cm/2.0;
           G4double startAngleOfTheHolder = 0*deg;
           G4double spanningAngleOfTheHolder = 360*deg;

           // Define volumes.
           G4Tubs* sourceHolderLeadCo57 = new G4Tubs("sourceHolderLeadCo57",
               innerRadiusOfTheLeadCo57,
               outerRadiusOfTheLeadCo57,
               heightOfTheLeadCo57,
               startAngleOfTheHolder,
               spanningAngleOfTheHolder);

           // Position of the extra lead for Co 57.
           G4ThreeVector sourceHolderCo57LeadPos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator+heightOfTheLeadInside1+2.0*heightOfTheLeadCo57);

           // Define logical and physical volume.
           sourceHolderCo57LeadLogic = new G4LogicalVolume(sourceHolderLeadCo57, Pb, "sourceHolderLeadCo57",0,0,0);
           sourceHolderCo57LeadPhys = new G4PVPlacement(0, sourceHolderCo57LeadPos, sourceHolderCo57LeadLogic,
               "sourceHolderCo57LeadPhys", experimentalHall_log, false, 0);
           
           // Is volume visible?
           G4VisAttributes* visiblesourceHolderCo57Lead = new G4VisAttributes(G4Colour(144,144,144)); 
           sourceHolderCo57LeadLogic->SetVisAttributes(visiblesourceHolderCo57Lead);
           // SetVisAttributes(G4VisAttributes::Invisible);//
           
           // Position of the lead cover if the source holder is Co 57.
           G4ThreeVector sourceHolderCoverPos = G4ThreeVector(holderAlX, holderAlY,
               holderAlZ+shiftCollimator+heightOfTheLeadInside1+heightOfTheSourceHolderCover2+heightOfTheHolderSecond/2.0);
         
           // Physical volume of the cover if the source holder is Co 57.
           sourceHolderCoverPhys = new G4PVPlacement(0, sourceHolderCoverPos, sourceHolderCoverLogic,
               "sourceHolderCoverPhys", experimentalHall_log, false, 0); 
         } else {
           
           // Position of the lead cover if the source holder is Am 241 or Cd 109.
           G4ThreeVector sourceHolderCoverPos = G4ThreeVector(holderAlX, holderAlY,
               holderAlZ+shiftCollimator+heightOfTheLeadInside1+heightOfTheSourceHolderCover2+heightOfTheHolderSecond/2.0-2.0*heightOfTheLeadInside3);

          // Physical volume of the cover id the source holder is Am 241 or Cd 109.
           sourceHolderCoverPhys = new G4PVPlacement(0, sourceHolderCoverPos, sourceHolderCoverLogic,
             "sourceHolderCoverPhys", experimentalHall_log, false, 0);
         }

         // Is volume visible?
         G4VisAttributes* sourceHolderCoverVolume = new G4VisAttributes(G4Colour(248,248,248)); 
         sourceHolderCoverLogic->SetVisAttributes(sourceHolderCoverVolume);
         // sourceHolderCoverLogic->SetVisAttributes(G4VisAttributes::Invisible);//
         // ------------------------------------------------------ //

         
         // ------------------------------------------------------ //
         // If collimator is disabled add extension (small cylindirical volume in front of the source holder).
         if (CGlobalParameters::enableCollimator == 0)
         {
           // Define dimensions.
           G4double innerRadiusOfTheHolderExtension = 0.4*cm;
           G4double outerRadiusOfTheHolderExtension = 0.75*cm;
           G4double heightOfTheHolderExtension = 1.0*cm/2.0;
           
           // Eliminate the shift in position due to the existence of the collimator.
           shiftCollimator = 2.0*heightOfTheHolderExtension;

           // Define volume
           G4Tubs* sourceHolderExtension = new G4Tubs("sourceHolderExtension",
               innerRadiusOfTheHolderExtension,
               outerRadiusOfTheHolderExtension,
               heightOfTheHolderExtension,
               startAngleOfTheHolder,
               spanningAngleOfTheHolder);
           
           // Position of the source holder extension.
           G4ThreeVector sourceHolderExtensionPos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator-heightOfTheHolderFirst-heightOfTheHolderExtension);
           
           // Define physical and logical volume.
           sourceHolderExtensionLogic = new G4LogicalVolume(sourceHolderExtension, Al, "sourceHolderExtension",0,0,0);
           sourceHolderExtensionPhys = new G4PVPlacement(0, sourceHolderExtensionPos, sourceHolderExtensionLogic,
               "sourceHolderExtensionPhys", experimentalHall_log, false, 0);
           
           // Is volume visible?
           G4VisAttributes* visiblesourceHolderExtension = new G4VisAttributes(G4Colour(248,248,248)); 
           sourceHolderExtensionLogic->SetVisAttributes(visiblesourceHolderExtension);
           //sourceHolderExtensionLogic->SetVisAttributes(G4VisAttributes::Invisible);
         }
         // ------------------------------------------------------ //
       }  
       break;
       
    case 2: // Source holder geometry is for Cs.
       { 
         G4cout << "Source holder for Cs." << G4endl;
         G4cout << " -----------------------------------------------------------" << G4endl;
         
         // ------------------------------------------------------ //
         // Define dimensions.
         G4double innerRadiusOfTheHolderCs137 = 0.2*cm;
         G4double outerRadiusOfTheHolderCs137 = 10.5*cm/2.0;
         G4double heightOfTheHolderCs137 = 8.4*cm/2.0;
         G4double startAngleOfTheHolder = 0*deg;
         G4double spanningAngleOfTheHolder = 360*deg;

         // Inner cut.
         G4double innerRadiusOfInnerCutCs137 = 0*cm;
         G4double outerRadiusOfInnerCutCs137 = 6.8*cm/2.0;
         G4double heightOfInnerCutCs137 = 7.3*cm/2.0;

         // Outer cut.
         G4double innerRadiusOfOuterCutCs137 = 8.0*cm/2.0;
         G4double outerRadiusOfOuterCutCs137 = 10.6*cm/2.0;
         G4double heightOfOuterCutCs137 = 7.3*cm/2.0;

         ReferenceFrame = heightOfTheHolderCs137+det_z/2.0;

         // Define coordinates for the position,
         G4double holderAlX = 0.0*cm;
         G4double holderAlY = 0.0*cm; 
         G4double holderAlZ = distanceToTheDetector+ReferenceFrame;

         // ShiftHolder. Needed below!.
         // 1.5*cm is the thickness of the collimator for Am, Co and Cd.
         // Shift depends on the existence of collimator.
         if ((CGlobalParameters::enableCollimator == 0) || (CGlobalParameters::enableCollimator == 2))
         {
           shiftCollimator = 1.0*cm;
         } else if (CGlobalParameters::enableCollimator == 1)
         {  
           shiftCollimator = 2.5*cm;
         } 
         
         // Define volumes.
         G4Tubs* InnerCut = new G4Tubs("InnerCut",
             innerRadiusOfInnerCutCs137,
             outerRadiusOfInnerCutCs137,
             heightOfInnerCutCs137,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* OuterCut = new G4Tubs("OuterCut",
             innerRadiusOfOuterCutCs137,
             outerRadiusOfOuterCutCs137,
             heightOfOuterCutCs137,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);

         G4Tubs* MainPieceCs = new G4Tubs("MainPieceCs",
             innerRadiusOfTheHolderCs137,
             outerRadiusOfTheHolderCs137,
             heightOfTheHolderCs137,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
 
         // Actually there is no rotation. It is needed for solid subtraction.
         G4RotationMatrix* zRot = new G4RotationMatrix;
         zRot->rotateZ(0*deg);
         zRot->rotateX(0*deg);
         zRot->rotateY(0*deg);

         // Subtract second lead from the first.
         G4ThreeVector zTransCs1(0,0,+(heightOfTheHolderCs137-heightOfInnerCutCs137));
         G4SubtractionSolid*  firstSubtractAlVolume = new G4SubtractionSolid("firstSubtractAlVolume", MainPieceCs, InnerCut, zRot, zTransCs1);

         // Subtract third lead from the remaining.
         G4ThreeVector zTransCs2(0,0,+(heightOfTheHolderCs137-heightOfInnerCutCs137));
         G4SubtractionSolid*  AlCs137Holder = new G4SubtractionSolid("AlCs137Holder", firstSubtractAlVolume, OuterCut, zRot, zTransCs2);
         
         // Position of the Cs 137 holder.
         G4ThreeVector AlCs137HolderPos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator);

         // Define physical and logical volume.
         AlCs137HolderLogic = new G4LogicalVolume(AlCs137Holder, Al, "AlCs137HolderLogic",0,0,0);        
         AlCs137HolderPhys = new G4PVPlacement(0, AlCs137HolderPos, AlCs137HolderLogic,
             "AlCs137HolderPhys", experimentalHall_log, false, 0);
         
         // Is volume visible?
         G4VisAttributes* visibleAlCs137Holder = new G4VisAttributes(G4Colour(248,248,248)); 
         AlCs137HolderLogic->SetVisAttributes(visibleAlCs137Holder);
         // ------------------------------------------------------ //
   
   
         // ------------------------------------------------------ //
         // If collimator is disabled add extension (small cylindirical volume in front of the source holder).
         if (CGlobalParameters::enableCollimator == 5)
         {
           // Define dimensions.
           G4double innerRadiusOfTheHolderExtensionCs137 = 0.3*cm;
           G4double outerRadiusOfTheHolderExtensionCs137 = 1.5*cm;
           G4double heightOfTheHolderExtensionCs137 = 1.0*cm/2.0;
           
           // Eliminate the shift in position due to the existence of the collimator.
           shiftCollimator = 2.0*heightOfTheHolderExtensionCs137;
           
           // Define volumes.
           G4Tubs* sourceHolderExtensionCs137 = new G4Tubs("sourceHolderExtensionCs137",
               innerRadiusOfTheHolderExtensionCs137,
               outerRadiusOfTheHolderExtensionCs137,
               heightOfTheHolderExtensionCs137,
               startAngleOfTheHolder,
               spanningAngleOfTheHolder);
           
           // Define positions.
           G4ThreeVector sourceHolderExtensionCs137Pos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator-heightOfTheHolderCs137-heightOfTheHolderExtensionCs137);
           
           // Define physical and logical volume.
           sourceHolderExtensionCs137Logic = new G4LogicalVolume(sourceHolderExtensionCs137, Al, "sourceHolderExtensionCs137",0,0,0);
           sourceHolderExtensionCs137Phys = new G4PVPlacement(0, sourceHolderExtensionCs137Pos, sourceHolderExtensionCs137Logic,
               "sourceHolderExtensionCs137Phys", experimentalHall_log, false, 0);
           
           // Is volume visible?
           G4VisAttributes* visiblesourceHolderExtensionCs137 = new G4VisAttributes(G4Colour(248,248,248)); 
           sourceHolderExtensionCs137Logic->SetVisAttributes(visiblesourceHolderExtensionCs137);
           // SetVisAttributes(G4VisAttributes::Invisible);//
         }
         // ------------------------------------------------------ //
      
         
         // ------------------------------------------------------ //
         // Lead inside the source holder Cs 137.
         // Define dimensions.
         G4double innerRadiusOfTheLeadInside1 = 0.0*mm;
         G4double outerRadiusOfTheLeadInside = 6.78*cm/2.0;
         G4double heightOfTheLeadInside1 = 6.0*cm/2.0;
         
         G4double outerRadiusOfTheLeadInside2 = 0.4*cm/2.0;
         G4double heightOfTheLeadInside2 = 2.5*cm/2.0;
         
         G4double outerRadiusOfTheLeadInside3 = 1.25*cm;
         G4double heightOfTheLeadInside3 = 0.5*cm/2.0;
         
         // Define volumes.
         G4Tubs* firstLeadInsideCs137 = new G4Tubs("firstLeadInsideCs137",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside,
             heightOfTheLeadInside1,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* secondLeadInsideCs137 = new G4Tubs("secondLeadInsideCs137",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside2,
             heightOfTheLeadInside2,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* thirdLeadInsideCs137 = new G4Tubs("thirdLeadInsideCs137",
             innerRadiusOfTheLeadInside1,
             outerRadiusOfTheLeadInside3,
             heightOfTheLeadInside3,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         // Subtract second lead from the first.
         G4ThreeVector zTransLeadCs1371(0,0,-(heightOfTheLeadInside1-heightOfTheLeadInside2)-0.01*mm);
         G4SubtractionSolid*  firstSubtractCs137Volume = new G4SubtractionSolid("firstSubtractCs137Volume", firstLeadInsideCs137, secondLeadInsideCs137, zRot, zTransLeadCs1371);
         
         // Subtract third lead from the remaining.
         G4ThreeVector zTransLeadCs1372(0,0,-heightOfTheLeadInside3);
         G4SubtractionSolid*  secondSubtractCs137Volume = new G4SubtractionSolid("secondSubtractCs137Volume", firstSubtractCs137Volume, thirdLeadInsideCs137, zRot, zTransLeadCs1372);
         
         // Position of the lead inside the source holder.
         G4ThreeVector LeadInsideCs137Pos = G4ThreeVector(holderAlX, holderAlY, holderAlZ+shiftCollimator);
         
         // Define physical and logical volumes.
         LeadInsideCs137Logic = new G4LogicalVolume(secondSubtractCs137Volume, Pb, "LeadInsideCs137Logic",0,0,0);
         LeadInsideCs137Phys = new G4PVPlacement(0, LeadInsideCs137Pos, LeadInsideCs137Logic,
             "LeadInsideCs137Phys", experimentalHall_log, false, 0);
         
         // Is volume visible?
         G4VisAttributes* visibleLeadInsideCs137 = new G4VisAttributes(G4Colour(1,0,0)); 
         LeadInsideCs137Logic->SetVisAttributes(visibleLeadInsideCs137);
         // ------------------------------------------------------ //
      
      
         // ------------------------------------------------------ //
         // Cap of the source holder for Cs 137.
         // Define dimensions.
         G4double innerRadiusOfTheSourceHolderCoverCs1371 = 6.8*cm/2.0;
         G4double outerRadiusOfTheSourceHolderCoverCs1371 = 8.02*cm/2.0;
         G4double heightOfTheSourceHolderCoverCs1371 = 2.5*cm/2.0;
         
         G4double innerRadiusOfTheSourceHolderCoverCs1372 = 0.0*mm;
         G4double outerRadiusOfTheSourceHolderCoverCs1372 = 8.0*cm/2.0;
         G4double heightOfTheSourceHolderCoverCs1372 = 3.8*cm/2.0;
         
         // Define volumes.
         G4Tubs* firstSourceHolderCs137Sub = new G4Tubs("firstSourceHolderCs137Sub",
             innerRadiusOfTheSourceHolderCoverCs1371,
             outerRadiusOfTheSourceHolderCoverCs1371,
             heightOfTheSourceHolderCoverCs1371,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         G4Tubs* firstSourceHolderCs137Cover = new G4Tubs("firstSourceHolderCs137Cover",
             innerRadiusOfTheSourceHolderCoverCs1372,
             outerRadiusOfTheSourceHolderCoverCs1372,
             heightOfTheSourceHolderCoverCs1372,
             startAngleOfTheHolder,
             spanningAngleOfTheHolder);
         
         // This is the shift needed to subtract correct part of the volume.        
         G4ThreeVector zTransCoverCs137(0,0,+heightOfTheSourceHolderCoverCs1371-heightOfTheSourceHolderCoverCs1372-0.01*mm);
         // This is the subtracted volume to be used ad a cap for the source holder.
         G4SubtractionSolid* sourceHolderCoverCs137 = new G4SubtractionSolid("sourceHolderCoverCs137", firstSourceHolderCs137Cover,
             firstSourceHolderCs137Sub, zRot, zTransCoverCs137);
      
         // Position of the cap.
         G4ThreeVector sourceHolderCoverCs137Pos = G4ThreeVector(holderAlX, holderAlY,
             holderAlZ+shiftCollimator+heightOfTheLeadInside1+heightOfTheSourceHolderCoverCs1372+heightOfTheLeadInside3/2.0);         
        
         // Define physical and logical volume.
         sourceHolderCoverCs137Logic = new G4LogicalVolume(sourceHolderCoverCs137, Al, "sourceHolderCoverCs137",0,0,0); 
         sourceHolderCoverCs137Phys = new G4PVPlacement(0, sourceHolderCoverCs137Pos, sourceHolderCoverCs137Logic,
             "sourceHolderCoverCs137Phys", experimentalHall_log, false, 0);
       
         // Is volume visible?
         G4VisAttributes* sourceHolderCoverCs137Volume = new G4VisAttributes(G4Colour(248,248,248)); 
         sourceHolderCoverCs137Logic->SetVisAttributes(sourceHolderCoverCs137Volume);
         // SetVisAttributes(G4VisAttributes::Invisible);//
         // ------------------------------------------------------ //
       }
       break;
  }
  // ******************************************************************** //
  // Source Holder ends here!
  // ******************************************************************** //


  // ******************************************************************** //
  // Sensitive detectors!
  // ******************************************************************** //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerSDname = "Ex/TrackerSD";
  ExTrackerSD* aTrackerSD = new ExTrackerSD(trackerSDname);
  SDman->AddNewDetector(aTrackerSD);
  CZTDet_log->SetSensitiveDetector(aTrackerSD);
  // ******************************************************************** //

  return experimentalHall_phys;
}
