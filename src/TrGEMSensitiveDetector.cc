
#include "TrGEMSensitiveDetector.hh"
#include "TrGEMAnalysis.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>

TrGEMSensitiveDetector::TrGEMSensitiveDetector(G4String name)
	:G4VSensitiveDetector(name), kickStart(false)
	{
  	collectionName.insert("TrGEMHitsCollection"); 	
	}

TrGEMSensitiveDetector::~TrGEMSensitiveDetector()
	{}

void TrGEMSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
	{
		hitsCollection = new TrGEMHitsCollection( SensitiveDetectorName, collectionName[0] );
		static int HCID = -1;
    
    if(HCID<0)
    	{ 
    		HCID = G4SDManager::GetSDMpointer()->GetCollectionID( collectionName[0] );
    	}
    HCE->AddHitsCollection( HCID, hitsCollection );
	}

G4bool TrGEMSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory* )
	{
		G4StepPoint *lastPoint = step->GetPostStepPoint();
		G4Track* track = step->GetTrack();
  		
   	G4int trackIndex=track->GetTrackID();
   	G4int parentIndex = track->GetParentID();
		
		G4double edep = step->GetTotalEnergyDeposit() ;
   	G4double nonionedep = step->GetNonIonizingEnergyDeposit() ;
   	edep -= nonionedep ;
   	
   	string particleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
   	
   	G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();
  	G4int chargePart = track->GetParticleDefinition()->GetPDGCharge();
  	G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
   	string volName = touchable->GetVolume(0)->GetName();
   	eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() ;
  	G4int stepNumber = track->GetCurrentStepNumber();
  	
		TrGEMHit *aHit = new TrGEMHit();
    aHit->SetTrackID(trackIndex);		
    aHit->SetCharge(chargePart);
    aHit->SetStepNumber(stepNumber);
    aHit->SetVolume(volName);
    aHit->SetKineticEnergy(track->GetDynamicParticle()->GetKineticEnergy());
    aHit->SetParentID(parentIndex);
    aHit->SetEdep(edep);
    aHit->SetPDG(pdg);
    aHit->SetCharge(chargePart);
    aHit->SetPos(lastPoint->GetPosition());
    aHit->SetMomentumDirection(lastPoint->GetMomentumDirection());  
    aHit->SetParticleName(particleName);  
    hitsCollection->insert(aHit);
    return true;	
	}
	
void TrGEMSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
	{}

