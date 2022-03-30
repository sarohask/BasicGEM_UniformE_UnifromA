
#include "TrGEMTrackingAction.hh"
#include "TrGEMAnalysis.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

TrGEMTrackingAction::TrGEMTrackingAction() 
	{;}

TrGEMTrackingAction::~TrGEMTrackingAction()
	{;}

void TrGEMTrackingAction::PreUserTrackingAction(const G4Track*) 
	{}

void TrGEMTrackingAction::PostUserTrackingAction(const G4Track*)
	{}
