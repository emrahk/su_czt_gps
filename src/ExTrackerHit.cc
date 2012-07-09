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
// $Id: ExTrackerHit.cc,v 1.9 2005/06/01 17:41:45 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<ExTrackerHit> ExTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExTrackerHit::ExTrackerHit()  : m_mom(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExTrackerHit::~ExTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExTrackerHit::ExTrackerHit(const ExTrackerHit& right)
  : G4VHit()
{
  m_trackID   = right.m_trackID;
  m_edep      = right.m_edep;
  m_totEdep   = right.m_totEdep;
  m_mom       = right.m_mom;
  m_parPrePEn = right.m_parPrePEn;
  m_parPostPEn = right.m_parPostPEn;

  m_pos       = right.m_pos;
  m_prePos    = right.m_prePos;
  m_momPrePDir  = right.m_momPrePDir;
  m_momPostPDir = right.m_momPostPDir;
 
  m_processName = right.m_processName;
  m_particleName = right.m_particleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const ExTrackerHit& ExTrackerHit::operator=(const ExTrackerHit& right)
{
  m_trackID   = right.m_trackID;
  m_edep      = right.m_edep;
  m_totEdep   = right.m_totEdep;
  m_mom       = right.m_mom;
  m_parPrePEn = right.m_parPrePEn;
  m_parPostPEn = right.m_parPostPEn;

  m_pos       = right.m_pos;
  m_prePos    = right.m_prePos;
  m_momPrePDir  = right.m_momPrePDir;
  m_momPostPDir = right.m_momPostPDir;
  
  m_processName = right.m_processName;
  m_particleName = right.m_particleName;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ExTrackerHit::operator==(const ExTrackerHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(m_pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExTrackerHit::Print()
{
  G4cout << "--> trackID: " << m_trackID 
	     << " Info: "
	     << " -->> Par:" << m_particleName
		 << " || -->> Pr=" << G4BestUnit(m_parPrePEn,"Energy")
		 << " || -->> Po=" << G4BestUnit(m_parPostPEn,"Energy")
		 << " || -->> DE=" << G4BestUnit(m_edep,"Energy")
		 << " || -->> TE=" << G4BestUnit(m_totEdep,"Energy")
		 << " || -->>  " << m_processName
		 << " || -->> Pr(" << G4BestUnit(m_prePos,"Length") << ")"
	     << " || -->> Po(" << G4BestUnit(m_pos,"Length") << ")"
		 << " || -->> Dr" << m_momPrePDir 
		 << " || -->> Do" << m_momPostPDir  
		 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

