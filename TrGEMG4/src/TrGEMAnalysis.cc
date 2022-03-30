
#include "GEMDPG/TrGEMG4/interface/TrGEMHit.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMAnalysis.hh"

#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4String.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>

#include <TVector3.h>

TrGEMAnalysis* TrGEMAnalysis::singleton = 0 ;

TrGEMAnalysis* TrGEMAnalysis::GetInstance() 
	{
		if ( singleton == 0 ) 
			{
    		static TrGEMAnalysis trGEMAnalysis;
    		singleton = &trGEMAnalysis;
  		}
  	return singleton;
	}

TrGEMAnalysis::~TrGEMAnalysis() 
	{}

TrGEMAnalysis::TrGEMAnalysis() 
	{
		file = 0;
	}
	
void TrGEMAnalysis::SetFileName(const G4String name)
	{
  	fileName = name;
	}

void TrGEMAnalysis::PrepareNewRun()
	{
		eventCounter = 0 ;
			
		fileName = fileName + ".root";
  	file = TFile::Open(fileName.c_str(),"RECREATE");
  	if(file) 
    	G4cout << "ROOT file " << fileName << " is created " << G4endl;
  	else 
  		{
    		G4cerr << "ROOT file has not been created!" << G4endl;
  		}
  	
  	// create tree and branch(es)
  	tHits = new TTree("HitsSteps", "Tree keeps output from G4Step Geant4 Simulation");
		tHits->Branch("EventID", &eventID_, "eventID_/I");
		tHits->Branch("TrackID", &trackID_, "trackID_/I");
   	tHits->Branch("Charge", &charge_);
   	tHits->Branch("StepNumber", &stepNum_, "stepNum_/I");
   	tHits->Branch("VolumeName", &volName_);
   	tHits->Branch("ParticleName", &particName_);
   	tHits->Branch("Edep", &eDep_);
   	tHits->Branch("PositionX", &positionX_);
   	tHits->Branch("PositionY", &positionY_);
   	tHits->Branch("PositionZ", &positionZ_);
   	tHits->Branch("PrimaryEnergy",  &primaryEnergy);  
   	tHits->Branch("PrimaryAngle",  &primaryAngle);
   	tHits->Branch("PrimaryPositionX",  &primaryPosX);
   	tHits->Branch("PrimaryPositionY",  &primaryPosY);
   	tHits->Branch("PrimaryPositionZ",  &primaryPosZ);
	}

void TrGEMAnalysis::PrepareNewEvent()
	{
  	//Reset variables relative to this event
  	eventID_ = 0;
   	stepNum_ = -1.0 ;
 	}

void TrGEMAnalysis::EndOfEvent()
	{
  	tHits->Fill() ;
	}

void TrGEMAnalysis::FillPrimary()
	{
  	tEvent->Fill();
	}

void TrGEMAnalysis::EndOfRun()
	{
  	// Writing and closing the ROOT file
  	G4cout << "ROOT: files writing..." << G4endl;
  	file->Write();
  	G4cout << "ROOT: files closing..." << G4endl;
  	file->Close();
  	G4cout << "ROOT: files closed..." << G4endl;
  	delete file ;
	}

void TrGEMAnalysis::SavePrimary(G4double aPrimEnergy, G4double aPrimAngle, G4int aPrimParticlePDG, G4double aPrimMomX, G4double aPrimMomY, G4double aPrimMomZ, G4double aPrimPositionX, G4double aPrimPositionY, G4double aPrimPositionZ)
	{
		primaryEnergy = aPrimEnergy / MeV;
  	primaryAngle = aPrimAngle ;
  	primaryPDG = aPrimParticlePDG ;
  	primaryMomX = aPrimMomX ;
  	primaryMomY = aPrimMomY ;
  	primaryMomZ = aPrimMomZ ;
  	primaryPosX = aPrimPositionX ;
  	primaryPosY = aPrimPositionY ;
  	primaryPosZ = aPrimPositionZ ;
	}
	
void TrGEMAnalysis::SetValuesForInteractingEvents( G4int evtNo, TrGEMHit* hit)
	{
    eventID_ = evtNo;
    trackID_ = hit->GetTrackID();
    charge_ = hit->GetCharge();
   	stepNum_ = hit->GetStepNumber();
   	volName_ = hit->GetVolume();
   	eDep_ = hit->GetEdep();
    pdgEncoding_ = hit->GetPDG();
    kineticEnergy_ = hit->GetKineticEnergy() ;
    positionX_ = hit->GetPos().getX() ;
    positionY_ = hit->GetPos().getY() ;
    positionZ_ = hit->GetPos().getZ() ;
    momDirectionX_ = hit->GetMomentumDirection().getX() ; 
   	momDirectionY_ = hit->GetMomentumDirection().getY() ; 
   	momDirectionZ_ = hit->GetMomentumDirection().getZ() ; 
   	momentumX_ = hit->GetMomentum().getX() ; 
   	momentumY_ = hit->GetMomentum().getY() ; 
   	momentumZ_ = hit->GetMomentum().getZ() ;
   	particName_ = hit->GetParticleName();	 	
	}

void TrGEMAnalysis::SetValuesForNonInteractingEvents(G4int evtNo)
	{
  	eventID_ = evtNo;
	}

void TrGEMAnalysis::SetKickStart(G4bool someKickStart) 
	{
		aKickStart = someKickStart ;
	}
	


	
