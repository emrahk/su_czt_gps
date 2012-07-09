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
// $Id: ExPrimaryGeneratorAction.hh,v 1.4 2002/01/09 17:23:48 ranjard Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#ifndef ExPrimaryGeneratorAction_h
#define ExPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"

class G4GeneralParticleSource;
class G4Event;

class ExPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExPrimaryGeneratorAction();
    ~ExPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:

    G4GeneralParticleSource* m_particleGun;
};

#endif

