
#include "GEMDPG/TrGEMG4/interface/TrGEMHit.hh"

TrGEMHit::TrGEMHit()
	: G4VHit(), fPDG(0), fStepNum(-1), fTrackID(-1), fEdep(0.0), fPos(0)   
	{}

TrGEMHit::~TrGEMHit()
	{}


