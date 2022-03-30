
#include "TrGEMStackingAction.hh"
#include "TrGEMAnalysis.hh"

#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

TrGEMStackingAction::TrGEMStackingAction()
	{}

TrGEMStackingAction::~TrGEMStackingAction() 
	{}

G4ClassificationOfNewTrack TrGEMStackingAction::ClassifyNewTrack( const G4Track*) 
	{
   	// always "urgent" in current applications
   	G4ClassificationOfNewTrack result( fUrgent );
   	return result;
	}

