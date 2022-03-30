
#include "GEMDPG/TrGEMG4/interface/TrGEMRunAction.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMAnalysis.hh"

#include "G4Run.hh"

TrGEMRunAction::TrGEMRunAction()
	{}

TrGEMRunAction::~TrGEMRunAction()
	{}

void TrGEMRunAction::BeginOfRunAction(const G4Run* aRun)
	{
   	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl ;
   	TrGEMAnalysis::GetInstance()->PrepareNewRun() ;
	}

void TrGEMRunAction::EndOfRunAction(const G4Run*)
	{
   	TrGEMAnalysis::GetInstance()->EndOfRun() ;
	}

