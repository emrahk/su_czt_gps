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
// $Id: ExTrackerHit.hh,v 1.7 2003/05/28 09:54:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExTrackerHit_h
#define ExTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExTrackerHit : public G4VHit
{
  public:

    ExTrackerHit();
    ~ExTrackerHit();
    ExTrackerHit(const ExTrackerHit&);
    const ExTrackerHit& operator=(const ExTrackerHit&);
    G4int operator==(const ExTrackerHit&) const;
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
    void Draw();
    void Print();
    
  public:
    
    void SetTrackID     (G4int track)       { m_trackID = track; };
    void SetEdep        (G4double de)       { m_edep = de; };
    void SetTotEdep     (G4double totde)    { m_totEdep = totde; }
    void SetMom         (G4ThreeVector mo)  { m_mom = mo;}
    void SetParPrePEn   (G4double en)       { m_parPrePEn = en; }
    void SetParPostPEn  (G4double en)       { m_parPostPEn = en; }
    
    void SetPos         (G4ThreeVector xyz) { m_pos = xyz; };
    void SetPrePos	  (G4ThreeVector xyz) { m_prePos = xyz; }
    void SetMomPrePDir  (G4ThreeVector xyz) { m_momPrePDir = xyz; }
    void SetMomPostPDir (G4ThreeVector xyz) { m_momPostPDir = xyz; }
    
    void SetProcessName (G4String pname)    { m_processName = pname; }
    void SetParticleName(G4String paname)   { m_particleName = paname; }
    
    void SetInitEnergy   (G4double energy)  { m_initEnergy = energy;}
    void SetTotDepEnergy (G4double energy)  { m_totDepEnergy = energy;}
    void SetEventType    (G4int eventtype)  { m_eventType = eventtype;}
    void SetGammaHitNumber (G4int hitnum)   { m_hitNumber = hitnum;}
    void SetGammaHitPattern (G4int pattern) { m_pattern = pattern;}
    void SetInitParID (G4int initparid)     { m_initParID = initparid; }
    
    void SetIfCrossed (G4bool crossed)	  { m_crossed = crossed; }
    void SetIfInteracted (G4bool interacted){ m_interacted = interacted; }
    
    G4int         GetTrackID()     { return m_trackID; };
    G4double      GetEdep()        { return m_edep; };  
    G4double      GetTotEdep()     { return m_totEdep; } 
    G4ThreeVector GetMoM()         { return m_mom; }
    G4double      GetParPrePEn()   { return m_parPrePEn; }
    G4double      GetParPostPEn()  { return m_parPostPEn; }
    
    G4ThreeVector GetPos()         { return m_pos; };
    G4ThreeVector GetPrePos()      { return m_prePos; }
    G4ThreeVector GetMomPrePDir()  { return m_momPrePDir; }
    G4ThreeVector GetMomPostPDir() { return m_momPostPDir; }
    
    G4String      GetProcessName() { return m_processName; }
    G4String      GetParticleName(){ return m_particleName; }
    
    G4double      GetInitEnergy()  { return m_initEnergy;}
    G4double		GetTotDepEnergy() { return m_totDepEnergy;}
    G4int		    GetEventType()     { return m_eventType;}
    G4int         GetGammaHitNumber() { return m_hitNumber;}
    G4int			GetGammaHitPattern() { return m_pattern;}
    G4int			GetInitParID() { return m_initParID; }
    
    G4bool	    IsCrossed() { return m_crossed; }
    G4bool	IsInteracted() { return m_interacted; }
    
  private:
    
    G4int         m_trackID;
    G4double      m_edep;
    G4ThreeVector m_mom;
    G4double      m_parPrePEn;
    G4double      m_parPostPEn;
    G4double      m_totEdep;
    
    G4ThreeVector m_pos;
    G4ThreeVector m_prePos;
    G4ThreeVector m_momPrePDir;
    G4ThreeVector m_momPostPDir;
    
    G4String      m_processName;
    G4String      m_particleName;
    
    G4double      m_initEnergy;
    G4double		m_totDepEnergy;
    G4int		    m_eventType;
    G4int         m_hitNumber;
    G4int			m_pattern;
    G4int			m_initParID;
    
    G4bool		m_interacted;
    G4bool		m_crossed;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<ExTrackerHit> ExTrackerHitsCollection;

extern G4Allocator<ExTrackerHit> ExTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
inline void* ExTrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExTrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
inline void ExTrackerHit::operator delete(void *aHit)
{
  ExTrackerHitAllocator.FreeSingle((ExTrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
