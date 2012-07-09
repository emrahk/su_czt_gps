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
// $Id: ExTrackerSD.hh,v 1.6 2002/01/09 17:24:09 ranjard Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExTrackerSD_h
#define ExTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExTrackerHit.hh"

#include "Randomize.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExTrackerSD : public G4VSensitiveDetector
{
  public:
    ExTrackerSD(G4String);
    ~ExTrackerSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
    
  private:
    ExTrackerHitsCollection* trackerCollection;
    
  private:
    G4double m_totEdep;
    G4int m_eventType;
    G4double m_firstAngel;
    G4double m_firstEdep;
    G4double m_firstElost;
    G4double m_allComptEn;
    
    G4bool m_bCrossed;
    G4bool m_bInteracted;
    
    // FILE *fpPhotoEM;
    // std::ofstream fpPhotoEM;
    
    //G4RandGauss *m_gaussRand;
    
    G4int m_nInitParTrackID;
    G4long m_nHitPattern;
    G4int m_nGammaHit;
    
    G4double m_fInitParticleEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

