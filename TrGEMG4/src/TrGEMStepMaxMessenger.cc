
#include "GEMDPG/TrGEMG4/interface/TrGEMStepMaxMessenger.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMStepMax.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

StepMaxMessenger::StepMaxMessenger(StepMax* stepM)
  : fStepMax(stepM)
	{ 
  	fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/testem/stepmax",this);
  	fStepMaxCmd->SetGuidance("Set max allowed step length");
  	fStepMaxCmd->SetParameterName("mxStep",false);
  	fStepMaxCmd->SetRange("mxStep>0.");
  	fStepMaxCmd->SetUnitCategory("Length");
  	fStepMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	}

StepMaxMessenger::~StepMaxMessenger()
	{
  	delete fStepMaxCmd;
	}

void StepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
	{ 
  	if (command == fStepMaxCmd)
    	{ 
    		fStepMax->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
    	}
	}

