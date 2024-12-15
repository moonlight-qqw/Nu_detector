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
// * work  make  any representation or  warranty, express or implied, *
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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Analysis.hh"
#include <math.h>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="opticalphoton");
  //fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  //fParticleGun->SetParticleEnergy(1.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//fParticleGun->SetParticlePosition(G4ThreeVector(-1500, 0, 0));
	
	// random distribution inside a spherical scintillator volume
	
	
	
	G4int option = 1; // 0 -- random distr    1 -- uniform distr 
	
	
	G4double RADIUS = 620; // radius of scintillator sphere
	if (option){
	G4double r = 600*((float) rand()/RAND_MAX); // 0 -- 630 (scint sphere rad) 
	G4double phi = 2*M_PI*((float) rand()/RAND_MAX); // 0 -- 2pi
	G4double theta = M_PI*((float) rand()/RAND_MAX) - M_PI/2; // -pi/2 -- pi/2
	
	G4double x0 = r*sin(theta)*cos(phi);
	G4double y0 = r*sin(theta)*sin(phi);
	G4double z0 = r*cos(theta);
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
	G4double dist = sqrt(x0*x0+y0*y0+z0*z0);
	
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-x0,-y0,-z0));
  	fParticleGun->SetParticleEnergy(1.*MeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);
	analysisManager->FillH1(2, dist);
	} else {
	 G4double phi = 2*M_PI*((float) rand()/RAND_MAX);
	 G4double costheta = 2*((float) rand()/RAND_MAX)-1;
	 G4double u = ((float) rand()/RAND_MAX);
	 G4double theta = acos(costheta);
	 G4double r = RADIUS * pow(u,1/3);
	 G4double x0 = r* sin(theta) * cos(phi);
         G4double y0 = r* sin(theta) * sin(phi);
         G4double z0 = r * cos(theta);
         G4double dist = sqrt(x0*x0+y0*y0+z0*z0);
         analysisManager->FillH1(2, dist);
         }
	    //phi = np.random.uniform(0, 2 * np.pi, size=PTS_COUNT - 1)
    //costheta = np.random.uniform(-1, 1, size=PTS_COUNT - 1)
    //u = np.random.random(size=PTS_COUNT - 1)
    //theta = np.arccos(costheta)
    //r = RADIUS * np.cbrt(u)
	// 
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
/*
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  G4double size = 0.8;
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = -0.5 * envSizeZ;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  

  fParticleGun->GeneratePrimaryVertex(anEvent);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


