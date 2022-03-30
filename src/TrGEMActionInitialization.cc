
#include "TrGEMActionInitialization.hh"
#include "TrGEMPrimaryGeneratorAction.hh"
#include "TrGEMRunAction.hh"
#include "TrGEMEventAction.hh"
#include "TrGEMSteppingAction.hh"
#include "TrGEMStackingAction.hh"
#include "TrGEMTrackingAction.hh"

#include "G4String.hh"

TrGEMActionInitialization::TrGEMActionInitialization(char* particleName)
	: G4VUserActionInitialization(), partName(particleName)
	{}

TrGEMActionInitialization::~TrGEMActionInitialization()
	{}

void TrGEMActionInitialization::BuildForMaster() const 
	{
  	TrGEMRunAction* runAction = new TrGEMRunAction;
		SetUserAction(runAction);
	}

void TrGEMActionInitialization::Build() const 
	{
  	SetUserAction(new TrGEMPrimaryGeneratorAction(partName));
	
  	TrGEMRunAction* runAction = new TrGEMRunAction;
  	SetUserAction(runAction);

  	TrGEMEventAction* eventAction = new TrGEMEventAction();
  	SetUserAction(eventAction);

  	TrGEMSteppingAction* steppingAction = new TrGEMSteppingAction();
  	SetUserAction(steppingAction);

  	TrGEMStackingAction* stackingAction = new TrGEMStackingAction();
  	SetUserAction(stackingAction);
  
  	TrGEMTrackingAction* trackingAction = new TrGEMTrackingAction();
  	SetUserAction(trackingAction);
	}
