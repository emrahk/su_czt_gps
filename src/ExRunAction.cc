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
// * work  make  any ExRunAction or  warranty, express or implied, *
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
//
// $Id: ExRunAction.cc,v 1.18 2006/06/29 17:49:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "ExGlobalParameters.hh"

#include "Randomize.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRunAction::ExRunAction()
{
  theSDName.push_back(G4String("Ex/TrackerSD"));
  m_nRunID = 0;
  
  timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
ExRunAction::~ExRunAction()
{
  theSDName.clear();
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
void ExRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  if ( m_nRunID == 0 ) G4cout << G4endl;
  G4cout << "========================== Run " << m_nRunID << " start. =============================" << G4endl;
  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  //initialize cumulative quantities
  //
  timer->Start();
  
  HepRandom::setTheSeed(time(0));
}

  
void ExRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  
  ExRun* reRun = (ExRun*)aRun;
  
  G4long simuNum = reRun->GetNumberOfEvent();
  G4long hittedNum = reRun->m_nHittedEventNum;
  G4long crossedNum = reRun->m_nCrossedEventNum;
  G4long peakNum = reRun->m_nPeakEventNum; 
  
  G4cout << "Passed Time: " << *timer << "|" << G4endl;
  G4cout << "Total simulated Events: "<< simuNum << "|"
    << "Events Passed the detector: " << crossedNum << "| "
    << "Events Interacted: " << hittedNum << "| "
    << "Photopeak Events: " << peakNum
    << G4endl;

  
  m_nRunID++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4Run* ExRunAction::GenerateRun()
{
  ExRun * userRun = new ExRun(theSDName);
  userRun->SetRunID(m_nRunID);
  return userRun;
}
