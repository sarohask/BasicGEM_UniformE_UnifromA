
#ifndef TrGEMHit_h
#define TrGEMHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "globals.hh"

using namespace std;

class TrGEMHit : public G4VHit
	{
    public:
    	TrGEMHit(); 
    	virtual ~TrGEMHit();
    	
    public:
    	inline void SetEdep(G4double de)  																	{ fEdep = de; }
    	inline G4double GetEdep() 																					{ return fEdep; }
    	inline void AddEdep(G4double dep) 																	{ fEventEdep += dep; }
    	inline G4double GetEventEdep() 																			{ return fEventEdep; }
    	
    	inline void SetPDG(G4int pdg) 																			{ fPDG = pdg; }
    	inline G4int GetPDG() 																							{ return fPDG; }
    	
    	inline void SetCharge(G4int aCharge) 																{ fCharge = aCharge; }
    	inline G4int GetCharge() 																						{ return fCharge; }
    	
    	inline void SetParentID(G4int parentID) 														{ fParentID = parentID; }
    	inline G4int GetParentID() 																					{ return fParentID; }
    	
    	inline void SetTrackID(G4int trkID) 																{ fTrackID = trkID; }
    	inline G4int GetTrackID() 																					{ return fTrackID; }
    	
    	inline void SetKineticEnergy(G4double kin) 													{ fKin = kin; }
    	inline G4double GetKineticEnergy() 																	{ return fKin; }
    	
    	inline void SetPos(G4ThreeVector xyz) 															{ fPos =xyz; }
    	inline G4ThreeVector GetPos() 																			{ return fPos; }
    	
    	inline void SetStepNumber(G4int stepNum) 														{ fStepNum = stepNum; }
    	inline G4int GetStepNumber() 																				{ return fStepNum; }	
    	
    	inline void SetInterLength(G4double interLen) 											{ fInterLength = interLen; }
    	inline G4double GetInterLength() 																		{ return fInterLength; }																		
    	
    	inline void SetMomentumDirection(G4ThreeVector moment) 							{ fMomDirection = moment; }
    	inline G4ThreeVector GetMomentumDirection() 												{ return fMomDirection; }
    	
    	inline void SetMomentum(G4ThreeVector momentum) 										{ fMomentum = momentum; }
    	inline G4ThreeVector GetMomentum() 																	{ return fMomentum; }
    	
    	inline void SetVolume(string vName) 																{ fVolName = vName; }
    	inline string GetVolume() 																					{ return fVolName; }
 
     	inline void SetProcessName(string procName) 	 											{ fProcName = procName; };
     	inline string GetProcessName()													 						{ return fProcName; }
     	
     	inline void SetParticleName(string particleName) 	 									{ fParticleName = particleName; };
     	inline string GetParticleName()													 						{ return fParticleName; }
     	
     	inline void SetNumInteraction (G4int numInter) 											{ fNumInteraction = numInter; }
    	inline G4int GetNumInteraction() 																		{ return fNumInteraction; }

		private:
    	G4int 															fPDG;
    	G4int 															fCharge;
    	G4int 															fStepNum;    	
    	G4int 															fParentID;
    	G4int 															fTrackID;
    	G4int 															fNumInteraction;
    	G4double														fInterLength ;
    	G4double 														fEdep;
    	G4double 														fEventEdep;
    	G4double 														fKin;
    	G4ThreeVector 											fPos; 
    	G4ThreeVector 											fMomDirection;
    	G4ThreeVector												fMomentum;
  
    	string 															fVolName;
 	 		string 															fProcName;
 	 		string 															fParticleName;
	};

typedef G4THitsCollection<TrGEMHit>TrGEMHitsCollection;


#endif

