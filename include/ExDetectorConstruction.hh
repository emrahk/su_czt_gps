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
// $Id: ExDetectorConstruction.hh,v 1.5 2002/01/09 17:23:48 ranjard Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#ifndef ExDetectorConstruction_H
#define ExDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4tgrMessenger;

#include "G4VUserDetectorConstruction.hh"

class ExDetectorConstruction : public G4VUserDetectorConstruction
{

  public:
    
    ExDetectorConstruction();
    ~ExDetectorConstruction();
    
    G4VPhysicalVolume* Construct();
    
  private:
    
    // Logical volumes
    // Experimental hall and the detector.
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* CZTDet_log;
    G4LogicalVolume* CZTDead_log;
    
    
    // Physical volumes
    // Experimental hall and the detector.
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* CZTDet_phys;
    G4VPhysicalVolume* CZTDead_phys;
    
  
    // G4AssemblyVolume;
    // This is used if there is an array of detectors.
    // Refer to EXDetectocConstruction.cc file.  
    G4tgrMessenger* messenger;
    
   
    // -----------------------------------------------------------//
    // Physical and logical volume for collimator.
    // Start!

    // Tungsten.
    G4LogicalVolume* collTungstenLogic;
    G4VPhysicalVolume* collTungstenPhys; 

    // Lead.
    G4LogicalVolume* collLeadLogic;
    G4VPhysicalVolume* collLeadPhys;

    //Aluminum Big.
    G4LogicalVolume* collAluminumBigLogic;
    G4VPhysicalVolume* collAluminumBigPhys;

    //Aluminum Small.
    G4LogicalVolume* collAluminumSmallLogic;
    G4VPhysicalVolume* collAluminumSmallPhys;

    //Aluminum Cover.
    G4LogicalVolume* collAluminumCoverLogic;
    G4VPhysicalVolume* collAluminumCoverPhys;
    // End!
    // -----------------------------------------------------------//


    // -----------------------------------------------------------//
    // Physical and logical volume for source holder. Am, Cd and Co.
    // Start!
 
    // Aluminum Source Holder.
    G4LogicalVolume* holderUnionVolumeLogic;
    G4VPhysicalVolume* holderUnionVolumePhys;

    // Aluminum Source Holder Extension if collimator is not present.
    G4LogicalVolume* sourceHolderExtensionLogic;
    G4VPhysicalVolume* sourceHolderExtensionPhys;

    // Lead pieces inside the source holder.
    G4LogicalVolume* firstSubtractLogic;
    G4VPhysicalVolume* firstSubtractPhys; 

    // Cover of the source holder.
    G4LogicalVolume* sourceHolderCoverLogic;
    G4VPhysicalVolume* sourceHolderCoverPhys;
    
    // Lead piece for Co 57 closer to the source holder cover.
    G4LogicalVolume* sourceHolderCo57LeadLogic;
    G4VPhysicalVolume* sourceHolderCo57LeadPhys;
    
    // Al box Cover.
    G4LogicalVolume* AlBoxCoverLogic;
    G4VPhysicalVolume* AlBoxCoverPhys;
    // End!
    // -----------------------------------------------------------//


    // -----------------------------------------------------------// 
    // Physical and logical volume for source holder. Cs.
    // Start!
    
    // Al box Cs 137 holder/
    G4LogicalVolume* AlCs137HolderLogic;
    G4VPhysicalVolume* AlCs137HolderPhys;

    // Al Source Holder Extension for Cs 137 if collimator is not present.
    G4LogicalVolume* sourceHolderExtensionCs137Logic;
    G4VPhysicalVolume* sourceHolderExtensionCs137Phys;

    //Lead inside the Source Holder for Cs 137 if collimator is not present.
    G4LogicalVolume* LeadInsideCs137Logic;
    G4VPhysicalVolume* LeadInsideCs137Phys;

    // Cover of the source holder for Cs 137.
    G4LogicalVolume* sourceHolderCoverCs137Logic;
    G4VPhysicalVolume* sourceHolderCoverCs137Phys;
    // End!
    // -----------------------------------------------------------//

    // Test collimator.
    G4LogicalVolume* collimatorTestLogic;
    G4VPhysicalVolume* collimatorTestPhys;

};
#endif

