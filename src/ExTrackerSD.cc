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
// $Id: ExTrackerSD.cc,v 1.7 2003/05/28 09:54:10 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

#include "ExGlobalParameters.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExTrackerSD::ExTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");

  //CreateDirectory("output",NULL);
  //fpPhotoEM=fopen("output\\thinCZT.gt4","wt");
  //fpPhotoEM.open("output\\thinCZT.gt4");

  //m_gaussRand = new G4RandGauss(HepRandom::getTheEngine());
  //G4RandGauss::shoot(0,0,sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExTrackerSD::~ExTrackerSD()
{
	//fclose(fpPhotoEM);
	//fpPhotoEM.close();
	//delete m_gaussRand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExTrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new ExTrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 

  m_totEdep = 0;
  m_eventType = 0;

  m_bCrossed = false;
  m_bInteracted = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ExTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4Track* track = aStep->GetTrack();             // Get my track	

  /*if ( track->GetDefinition()->GetParticleName() == 'e-' )
  {
	 track->SetTrackStatus(fStopAndKill);
	 return true;
  }*/

  G4double jm_edep = aStep->GetTotalEnergyDeposit();

  m_totEdep += jm_edep;
  m_bCrossed = true;
	
  //---------------------------------------------------------------
  // THIS IS THE POSITION TO GET ALL THE DATA IN EACH STEP
  // Information from G4Step
  G4ThreeVector jm_deltaposi = aStep->GetDeltaPosition();
  G4double      jm_deltatime = aStep->GetDeltaTime();
  G4ThreeVector jm_deltamome = aStep->GetDeltaMomentum();
  G4double      jm_deltaener = aStep->GetDeltaEnergy();
	
  // Information from G4StepPoint
  //   Prepoint
  G4StepPoint * jm_PrePoint = aStep->GetPreStepPoint();
  G4VPhysicalVolume * jm_PrePV  =  jm_PrePoint->GetPhysicalVolume();
  G4String      jm_PrePVname    =  jm_PrePV->GetName();
  G4ThreeVector jm_PrePPos      =  jm_PrePoint->GetPosition();
  G4ThreeVector jm_PrePMomDir   =  jm_PrePoint->GetMomentumDirection();
  G4Material *  jm_PrePMat      =  jm_PrePoint->GetMaterial();
  G4String      jm_PrePMatname  =  jm_PrePMat->GetName();
  G4String      jm_PrePProcname =  "null"; //track->GetCreatorProcess()->GetProcessName();
  G4double      jm_PrePTotEn    =  jm_PrePoint->GetTotalEnergy();
  G4ThreeVector jm_PrePMom      =  jm_PrePoint->GetMomentum();
  G4double      jm_PrePEn       =  jm_PrePoint->GetKineticEnergy();

  if(aStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL) 
        {jm_PrePProcname = aStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();}
	
  //   PostPoint
  G4StepPoint * jm_PostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume * jm_PostPV = jm_PostPoint->GetPhysicalVolume();
  G4String      jm_PostPVname =   jm_PostPV->GetName();
  G4ThreeVector jm_PostPPos =     jm_PostPoint->GetPosition();
  G4ThreeVector jm_PostPMomDir =  jm_PostPoint->GetMomentumDirection();
  G4Material *  jm_PostPMat    =  jm_PostPoint->GetMaterial();
  G4String      jm_PostPMatname = jm_PostPMat->GetName();
  G4String      jm_PostPProcname = "null";
  G4double      jm_PostPTotEn    =  jm_PostPoint->GetTotalEnergy();
  G4ThreeVector jm_PostPMom      =  jm_PostPoint->GetMomentum();
  G4double      jm_PostPEn       =  jm_PostPoint->GetKineticEnergy();

  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL) 
        { jm_PostPProcname = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();}

  // Get TrackID
  G4int trackID = track->GetTrackID();

  // store first gamma interaction information
  if (!m_eventType) 
  {
	  m_nInitParTrackID = trackID;
	  m_fInitParticleEnergy = jm_PrePEn;
	
  }

  if (jm_edep > 0) m_bInteracted = true;

  if (jm_PostPProcname=="LowEnConversion" || jm_PostPProcname=="conv") { m_eventType=4;}

  if (!m_eventType && (jm_PostPProcname=="LowEnCompton" || jm_PostPProcname=="compt") ) 
  { m_eventType=1; m_firstAngel=jm_PostPMomDir.angle(jm_PrePMomDir); m_firstEdep=jm_edep; m_firstElost=jm_PrePEn-jm_PostPEn; }
  else if (!m_eventType && (jm_PostPProcname=="LowEnPhotoElec" || jm_PostPProcname=="phot") ) 
  { m_eventType=2; m_firstAngel=-2; m_firstEdep=jm_edep; m_firstElost=jm_PrePEn-jm_PostPEn; m_firstAngel=jm_PostPMomDir.angle(jm_PrePMomDir); }
  else if (!m_eventType && jm_PostPProcname=="LowEnRayleigh")
  { m_eventType=3; m_firstAngel=jm_PostPMomDir.angle(jm_PrePMomDir); }
  else if (!m_eventType)
  { m_eventType= 99;} 

  // store hits information
  ExTrackerHit* newHit = new ExTrackerHit();

  newHit->SetTrackID(trackID);
  newHit->SetParticleName(track->GetDefinition()->GetParticleName());

  newHit->SetEdep(jm_edep);
  newHit->SetTotEdep(m_totEdep);

  newHit->SetParPrePEn(jm_PrePEn);
  newHit->SetParPostPEn(jm_PostPEn);
  newHit->SetPos(jm_PostPPos);
  newHit->SetPrePos(jm_PrePPos);

  newHit->SetMom(jm_PostPMom);
  newHit->SetMomPrePDir(jm_PrePMomDir);
  newHit->SetMomPostPDir(jm_PostPMomDir);
  
  newHit->SetProcessName(jm_PostPProcname);

  trackerCollection->insert( newHit );
  
  if(CGlobalParameters::displayTrack) newHit->Print();
  //newHit->Draw();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
	ExTrackerHit* newHit = new ExTrackerHit();

	newHit->SetInitEnergy(m_fInitParticleEnergy);
	newHit->SetTotDepEnergy(m_totEdep);
	newHit->SetEventType(m_eventType);
	newHit->SetGammaHitNumber(m_nGammaHit);
	newHit->SetGammaHitPattern(m_nHitPattern);
	newHit->SetInitParID(m_nInitParTrackID);
	newHit->SetIfCrossed(m_bCrossed);
	newHit->SetIfInteracted(m_bInteracted);

	trackerCollection->insert( newHit );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

