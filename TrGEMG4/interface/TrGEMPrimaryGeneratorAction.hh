
#ifndef TrGEMPrimaryGeneratorAction_h
#define TrGEMPrimaryGeneratorAction_h 1

#include "GEMDPG/TrGEMG4/interface/TrGEMAnalysis.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include "TF1.h"

using namespace std;

class G4Event;

class TrGEMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
	{
  	public:
    	TrGEMPrimaryGeneratorAction(char*);            
    	~TrGEMPrimaryGeneratorAction();

  	public:
    	virtual void GeneratePrimaries(G4Event*);
    	
  	private:
    	char*           					partName;
    	G4GeneralParticleSource* 	fParticleGun; 
    	G4double									primaryEne;
    	G4int           					primPDG;
    	G4double									primEnergy ;
    	G4double									primAngle ;
    	G4double									primPosX ;
    	G4double									primPosY ;
    	G4double									primPosZ ;
    	G4double									primMomX ;
    	G4double									primMomY ;
    	G4double									primMomZ ;
    	G4double									eneRange;
    	G4double									angRange;
    	G4double									angle;
	};

#endif
