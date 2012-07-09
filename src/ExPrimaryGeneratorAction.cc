
#include "ExPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "globals.hh"

#include "ExGlobalParameters.hh"
#include "Randomize.hh"

ExPrimaryGeneratorAction::ExPrimaryGeneratorAction()
{
   m_particleGun = new G4GeneralParticleSource();
  
  // Source position now will change according to the position of the collimator and/or source holder.
  // Added on 08/12/10 by Ozge.
  // defaultDistance changes according to the existence of the collimator and/or source holder.
  
  G4double sourcepos_x = 0.0*cm;
  G4double sourcepos_y = 0.0*cm;
  
  G4double defaultDistance; 

  // No collimator. No source holder.
  if ((!CGlobalParameters::enableCollimator) && (CGlobalParameters::enableSourceHolder) == 0)
  {
    defaultDistance = 5.0*cm;
  }
  // Collimator exist. Source holder is Am, Cd or Co.
  if ((CGlobalParameters::enableCollimator) && (CGlobalParameters::enableSourceHolder) == 1)
  {
    defaultDistance = 5.25*cm;
  }
  // No collimator. Source holder is Am, Cd or Co.
  if (!(CGlobalParameters::enableCollimator) && (CGlobalParameters::enableSourceHolder) == 1)
  {
    defaultDistance = 3.75*cm;
  }
  // Collimator exist. Source holder is Cs.  
   if ((CGlobalParameters::enableCollimator) && (CGlobalParameters::enableSourceHolder) == 2)
  {
    defaultDistance = 6.75*cm;
  }
  // No collimator. Source holder is Cs.     
  if (!(CGlobalParameters::enableCollimator) && (CGlobalParameters::enableSourceHolder) == 2)
  {
    defaultDistance = 5.25*cm;
  }

   G4double sourcepos_z = CGlobalParameters::sourceposZ + CGlobalParameters::distancetothedetector + CGlobalParameters::alboxcoverdistance + defaultDistance;
   G4SPSPosDistribution *posDist = m_particleGun->GetCurrentSource()->GetPosDist() ;
   posDist->SetCentreCoords(G4ThreeVector(sourcepos_x,sourcepos_y,sourcepos_z));

   G4cout << "sourcepos_x=" << sourcepos_x << "sourcepos_y=" << sourcepos_y << "sourcepos_z= " << sourcepos_z <<"mm"<< G4endl;

}

ExPrimaryGeneratorAction::~ExPrimaryGeneratorAction()
{
  delete m_particleGun;
}

void ExPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  m_particleGun->GeneratePrimaryVertex(anEvent) ;
}

