#ifndef TrGEMAnalysis_h
#define TrGEMAnalysis_h 1

#include "TrGEMHit.hh"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1.h"
#include "snprintf.h"

#include <iostream>
#include <fstream>

#include "globals.hh"
#include <vector>
#include <map>

using namespace std; 

class G4Run;
class G4Event;
class TFile;
class TTree;
class TString;

class TrGEMAnalysis 
	{
		public:
    	static TrGEMAnalysis* GetInstance();
    	~TrGEMAnalysis();
    	
    	void SetFileName(const G4String name);
    	void PrepareNewEvent();
    	void EndOfEvent();
    	void PrepareNewRun();
    	void EndOfRun();
    	void FillPrimary();
    	void SavePrimary(G4double aPrimEnergy, G4double aPrimAngle, G4int aPrimParticlePDG, G4double aPrimMomX, G4double aPrimMomY, G4double aPrimMomZ, G4double aPrimPositionX, G4double aPrimPositionY, G4double aPrimPositionZ);
    	void SetValuesForInteractingEvents ( G4int evtNo, TrGEMHit* hit);
    	void SetValuesForNonInteractingEvents(G4int evtNo);
    	void SetKickStart(G4bool aKickstart) ;
    	bool GetKickStart()																										{ return aKickStart ;}
    	
    	
  	private:
  		TrGEMAnalysis();
			static TrGEMAnalysis* singleton;
		
		public:
    	TFile*    				file;
    	TTree*    				tHits;
    	TTree*    				tEvent;
    	G4bool 						aKickStart ; 
    	G4int 						eventCounter ;
    	G4String 					fileName;	
    	G4int 						eventID_;
    	G4int 						trackID_;
  		G4int							charge_;	
   		G4int							stepNum_;
   		G4int							pdgEncoding_;
   		string		 				volName_;
   		string		 				particName_;
   		G4double 					eDep_;
  		G4double 					kineticEnergy_;
  		G4double 					positionX_ ;
  		G4double 					positionY_ ;
  		G4double 					positionZ_ ;
  		G4double					momDirectionX_ ; 
   		G4double					momDirectionY_ ; 
   		G4double					momDirectionZ_ ; 
   		G4double					momentumX_ ; 
   		G4double					momentumY_ ; 
   		G4double					momentumZ_ ;
   		
  		//tEvents  		
  		G4double 					primaryEnergy;
  		G4double 					primaryAngle;
  		G4double 					primaryPDG;
  		G4double 					primaryMomX;
  		G4double 					primaryMomY;
  		G4double 					primaryMomZ;
  		G4double 					primaryPosX;
  		G4double 					primaryPosY;
  		G4double 					primaryPosZ;
	};

#endif 
