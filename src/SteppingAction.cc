#include "SteppingAction.hh"
#include "ExGlobalParameters.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include <fstream>

SteppingAction::SteppingAction()
{;}

SteppingAction::~SteppingAction()
{;}

void SteppingAction::UserSteppingAction(const G4Step * theStep)
{
  G4Track * theTrack = theStep->GetTrack();
  G4StepPoint * thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume();
  G4String thePrePVname = thePrePV->GetName();
  
  G4ParticleDefinition * particleType = theTrack->GetDefinition();

  if(CGlobalParameters::stepview){
    std::ofstream ofs;
    ofs.open("output/stepping.txt",std::ios::app);
    if(thePrePVname(0,27) == "av_1_impr_1_CZTDet_log_pv_0" || thePrePVname(0,11) == "CZTDead_log"){
      ofs << theTrack->GetTrackID() << "  " << theTrack->GetParentID() << "  " << theTrack->GetGlobalTime() << "  " << theTrack->GetVertexPosition() << "  " << theTrack->GetVertexKineticEnergy()*1000 << "  " << theTrack->GetKineticEnergy()*1000 << "  " << (theTrack->GetVertexKineticEnergy()-theTrack->GetKineticEnergy())*1000 << "  ";
      ofs << theTrack->GetPosition() << " " << theTrack->GetMomentum() << "  ";
      if(particleType == G4Gamma::GammaDefinition()){ofs << "gamma" << G4endl;}
      if(particleType == G4Electron::ElectronDefinition()){ofs << "electron" << G4endl;}
    }
    ofs.close();
  }
  // then suspend the track
  //theTrack->SetTrackStatus(fSuspend);
  //return;

}
