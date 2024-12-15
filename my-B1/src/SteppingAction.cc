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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "Randomize.hh"
#include "Analysis.hh"

#include <fstream>
#include <cstdlib> 

using namespace std;
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //ofstream f;
  //int n_photons = 0;
  /*
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }
  */
  
  const G4String currentPhysicalName 
    = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  
  const G4String postPhysicalName 
    = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();  
    
  const G4String prePhysicalName 
    = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();  
    
  G4int CurrentTrackID;
  const G4String particleName
	= step->GetTrack()->GetDefinition()->GetParticleName();

  bool FirstPhoton = false;
  // get volume of the current step
  /*G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();*/

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;
  /*
  G4ThreeVector vec = step->GetPreStepPoint()->GetPosition();
  G4double rCurrent;
  G4double r = sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z());
  if (particleName == "opticalphoton" && !FirstPhoton){
  	//G4ThreeVector vec = step->GetPreStepPoint()->GetPosition();
  	//G4double r = sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z());
  	rCurrent = r;
  	CurrentTrackID = step -> GetTrack() -> GetTrackID();
  	G4cout << "Track: " << step -> GetTrack() -> GetTrackID() << G4endl;
  	G4cout << "RAD " << r << endl;
  	FirstPhoton = true;
  } else if (particleName == "opticalphoton" && FirstPhoton){
  	if (CurrentTrackID != step -> GetTrack() -> GetTrackID()){
  		CurrentTrackID = step -> GetTrack() -> GetTrackID();
  	  	//G4ThreeVector vec = step->GetPreStepPoint()->GetPosition();
  		//G4double r = sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z());
	  	G4cout << "Track: " << step -> GetTrack() -> GetTrackID() << G4endl;
	  	G4cout << "RAD " << r << endl;
	  	//analysisManager->FillH1(2, rCurrent);
  		FirstPhoton = false;
  		rCurrent = r;
  	} else { 
  		analysisManager->FillH1(2, rCurrent);
  		FirstPhoton = false;
  	*/	
  	
  	
  if (prePhysicalName == "pmt_glass"){// && particleName == "opticalphoton"){
  //G4cout << "/////////////////////" << G4endl;
  //G4cout << particleName << G4endl;
  
  fEventAction->nAbsPhotons++;
  G4double EdepStep = step->GetPreStepPoint()->GetKineticEnergy();
  
	  
 G4double prob_photon = 0.28;
 G4double rnd = (float) rand()/RAND_MAX;
  //G4double rnd = G4UniformRand();
  if (rnd < 0.28){
  	analysisManager->FillH1(0, EdepStep*1e06);
  	//G4cout << rnd << " " << "Energy: " << EdepStep <<  G4endl;
  	//G4cout << "/////////////////////" << G4endl;
  }
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  //G4cout << "Energy: " << EdepStep <<  G4endl;
  //analysisManager->FillH1(0, EdepStep*1e06);
  //G4cout << "/////////////////////" << G4endl;
  }
  
  //G4cout << "n_photons:  " << n_photons << G4endl;
  //f.open("output.txt", ios::app);
  //f << n_photons << endl;
  /*
  if (currentPhysicalName == "pmt"){
	  if (particleName == "opticalphoton"){
	  	G4cout << "//// photon ////" << G4endl;
	  	n_photons++;
		feventAction->nAbsPhotons++;
	    	//aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	  }
  }
  
  */
  
	//G4cout << "//////////////////////////////////" << G4endl;
	//G4cout << n_photons << G4endl;
	//G4cout << "//////////////////////////////////" << G4endl;
  // collect energy deposited in this step
  //G4double edepStep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
