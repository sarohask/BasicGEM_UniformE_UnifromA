
#include "TrGEMPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include <TMath.h>
#include <fstream>
#include "TRandom3.h"
#include "TH2.h"

TrGEMPrimaryGeneratorAction::TrGEMPrimaryGeneratorAction(char* partName_) 
	: partName(partName_)
	{
  	G4ParticleDefinition* particlen = G4ParticleTable::GetParticleTable()-> FindParticle(partName);
		
  	fParticleGun = new G4GeneralParticleSource();
  	fParticleGun->GetCurrentSource()->SetParticleDefinition(particlen);
  	
  	std::string temp = partName;
  	if (temp == "neutron") eneRange = 14;
  	else if (temp == "gamma") eneRange = 6;
  	else if (temp == "e-") eneRange = 5;
  	else if (temp == "e+") eneRange = 5;
  	else if (temp == "alpha") eneRange = 12;
  	else if (temp == "mu-") eneRange = 12;
  	else throw;
  	angRange = -90. ;
	}

TrGEMPrimaryGeneratorAction::~TrGEMPrimaryGeneratorAction()
	{
  	delete fParticleGun;
	}

void TrGEMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
	{
  	G4double xxx, zzz, gX(0.), gY(0.), gZ(0.);
  	G4double cosx(0.), cosy(0.), cosz(0.);
  	G4double tt, kk;
  
  	//set position  
  	gY = (2*G4UniformRand()-1)*260.;
  	xxx = G4UniformRand();
  	if (xxx <= 0.5)      
  		{ 
  			zzz = G4UniformRand();
  			gX = -zzz*((gY + 294.768589469)/11.2637724419) ;
  		}
  	if (xxx > 0.5)       
  		{ 
  			zzz = G4UniformRand();
  			gX = zzz*((gY + 294.768589469)/11.2637724419) ;
  		}
  	gZ = -18.57 ;	
  	
  	angle = (G4UniformRand()-1)*angRange;																						// - for even    + for odd
  		
  	if (gZ < 0.)
  		{
  	 		cosz = cos(angle * M_PI/180.0);	
  	 		
  	 		tt = G4UniformRand();
				if (tt < 0.5) 
					{
						cosy = abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;													//cos (beta)
					}
				if (tt >= 0.5)
					{
						cosy = -abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;												//cos (beta)
					}
				kk = G4UniformRand();
				if (kk < 0.5) 
					{
						cosx = abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
				if (kk >= 0.5)
					{
						cosx = -abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
			}
			
  	if (gZ > 0.)
  		{
  	 		cosz = -cos(angle * M_PI/180.0) ;
  	 		
  	 		tt = G4UniformRand();
				if (tt < 0.5) 
					{
						cosy = abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;													//cos (beta)
					}
				if (tt >= 0.5)
					{
						cosy = -abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;												//cos (beta)
					}
				kk = G4UniformRand();
				if (kk < 0.5) 
					{
						cosx = abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
				if (kk >= 0.5)
					{
						cosx = -abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}	
			}
			
		fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(gX*mm, gY*mm, gZ*mm));
		fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(cosx, cosy, cosz));
		
		primaryEne = TMath::Power(10, (G4UniformRand()-1)*eneRange);
  	fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(primaryEne*GeV);
 
  	//create vertex   
  	fParticleGun->GeneratePrimaryVertex(anEvent);
  	
  	primPDG = fParticleGun->GetParticleDefinition()->GetPDGEncoding() ;
  	primEnergy = fParticleGun->GetParticleEnergy() ;
  	if (gZ < 0.) { primAngle = acos(fParticleGun->GetParticleMomentumDirection().z())*180./M_PI;	}
  	if (gZ > 0.) { primAngle = acos((-1.)*fParticleGun->GetParticleMomentumDirection().z())*180./M_PI; }								//odd	
  	primMomX = fParticleGun->GetParticleMomentumDirection().x();
  	primMomY = fParticleGun->GetParticleMomentumDirection().y();
  	primMomZ = fParticleGun->GetParticleMomentumDirection().z();
  	primPosX = fParticleGun->GetParticlePosition().x() ;
  	primPosY = fParticleGun->GetParticlePosition().y() ;
  	primPosZ = fParticleGun->GetParticlePosition().z() ;
  	TrGEMAnalysis::GetInstance()->SavePrimary(primEnergy, primAngle, primPDG, primMomX, primMomY, primMomZ, primPosX, primPosY, primPosZ);
 }

