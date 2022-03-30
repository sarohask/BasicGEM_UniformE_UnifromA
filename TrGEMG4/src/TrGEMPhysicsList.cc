
#include "GEMDPG/TrGEMG4/interface/TrGEMPhysicsList.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMStepMax.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

TrGEMPhysicsList::TrGEMPhysicsList() : G4VModularPhysicsList()
	{
  	fConfig = G4LossTableManager::Instance()->EmConfigurator();
  	G4LossTableManager::Instance()->SetVerbose(1);
  	defaultCutValue = 1.*mm;
  	fCutForGamma     = defaultCutValue;
  	fCutForElectron  = defaultCutValue;
  	fCutForPositron  = defaultCutValue;
  	fCutForProton    = defaultCutValue;

  	fStepMaxProcess = new StepMax();

  	// Decay Physics is always defined
  	fDecayTrGEMPhysicsList = new G4DecayPhysics();

  	// EM physics
  	fEmName = G4String("emstandard");
  	fEmTrGEMPhysicsList = new G4EmStandardPhysics(1);
		
		SetVerboseLevel(1);
	}

TrGEMPhysicsList::~TrGEMPhysicsList()
	{
  	delete fDecayTrGEMPhysicsList;
  	delete fEmTrGEMPhysicsList;
  	for(size_t i=0; i<fHadronPhys.size(); ++i) 
  		{ 
  			delete fHadronPhys[i]; 
  		}
  	delete fStepMaxProcess;
	}

void TrGEMPhysicsList::ConstructParticle()
	{
  	fDecayTrGEMPhysicsList->ConstructParticle();
	}

void TrGEMPhysicsList::ConstructProcess()
	{
  	AddTransportation();
  	fEmTrGEMPhysicsList->ConstructProcess();
  	fDecayTrGEMPhysicsList->ConstructProcess();
  	for(size_t i=0; i<fHadronPhys.size(); ++i) 
  		{ 
  			fHadronPhys[i]->ConstructProcess(); 
  		}
  	AddStepMax();
	}

void TrGEMPhysicsList::AddTrGEMPhysicsList(const G4String& name)
	{
  	if (verboseLevel>1) 
  		{
    		G4cout << "TrGEMPhysicsList::AddTrGEMPhysicsList: <" << name << ">" << G4endl;
  		}

  	if (name == fEmName) 				{	return;} 
  	
  	else if (name == "emstandard_opt1") 
  		{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmStandardPhysics_option1();
  		} 
  	else if (name == "emstandard_opt2") 
  		{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmStandardPhysics_option2();
  		} 
  	else if (name == "emstandard_opt3") 
  		{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmStandardPhysics_option3();
  		} 
  	else if (name == "emstandard_opt4") 
  		{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmStandardPhysics_option4();
			} 
		else if (name == "emlivermore") 
			{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmLivermorePhysics();
  		} 
  	else if (name == "empenelope") 
  		{
    		fEmName = name;
    		delete fEmTrGEMPhysicsList;
    		fEmTrGEMPhysicsList = new G4EmPenelopePhysics();
			} 
		else if (name == "pai") 
			{
    		fEmName = name;
    		AddPAIModel(name);
  		} 
  	else if (name == "pai_photon") 
  		{ 
    		fEmName = name;
    		//AddPAIModel(name);
  		} 
  	else 
  		{
    		G4cout 	<< "TrGEMPhysicsList::AddTrGEMPhysicsList: <" << name << ">"
           			<< " is not defined"
           			<< G4endl;
  		}
	}

void TrGEMPhysicsList::AddStepMax()
	{
		fStepMaxProcess = new StepMax();
		
  	// Step limitation seen as a process
  	auto theParticleIterator = GetParticleIterator();
  	theParticleIterator->reset();
  	while ((*theParticleIterator)())
  		{
    		G4ParticleDefinition* particle = theParticleIterator->value();
    		G4ProcessManager* pmanager = particle->GetProcessManager();

    		if (fStepMaxProcess->IsApplicable(*particle))
    			{
      			pmanager ->AddDiscreteProcess(fStepMaxProcess);
    			}
  		}
	}

void TrGEMPhysicsList::SetCuts()
	{
  	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100.*eV,1e5);
  	if ( verboseLevel > 0 )
  		{
    		G4cout << "TrGEMPhysicsList::SetCuts:";
    		G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  		}

  	// set cut values for gamma at first and for e- second and next for e+,
  	// because some processes for e+/e- need cut values for gamma
  	SetCutValue(fCutForGamma, "gamma");
  	SetCutValue(fCutForElectron, "e-");
  	SetCutValue(fCutForPositron, "e+");
  	SetCutValue(fCutForProton, "proton");
  	if ( verboseLevel > 0 ) 
  		{ 
  			DumpCutValuesTable(); 
  		}
	}

void TrGEMPhysicsList::SetCutForGamma(G4double cut)
	{
  	fCutForGamma = cut;
  	SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
	}

void TrGEMPhysicsList::SetCutForElectron(G4double cut)
	{
  	fCutForElectron = cut;
  	SetParticleCuts(fCutForElectron, G4Electron::Electron());
	}

void TrGEMPhysicsList::SetCutForPositron(G4double cut)
	{
  	fCutForPositron = cut;
  	SetParticleCuts(fCutForPositron, G4Positron::Positron());
	}

void TrGEMPhysicsList::SetCutForProton(G4double cut)
	{
  	fCutForPositron = cut;
  	SetParticleCuts(fCutForProton, G4Proton::Proton());
	}

void TrGEMPhysicsList::AddPAIModel(const G4String& modname)
	{
		auto theParticleIterator = GetParticleIterator();
  	theParticleIterator->reset();
  	while ((*theParticleIterator)())
  		{
   	 		G4ParticleDefinition* particle = theParticleIterator->value();
    		G4String partname = particle->GetParticleName();
    		if(partname == "e-" || partname == "e+") 
    			{
      			NewPAIModel(particle, modname, "eIoni");
					} 
				else if(partname == "mu-" || partname == "mu+") 
					{
      			NewPAIModel(particle, modname, "muIoni");
					} 
				else if(partname == "proton" || partname == "pi+" || partname == "pi-") 
					{
      			NewPAIModel(particle, modname, "hIoni");
    			}
 			}
	}

void TrGEMPhysicsList::NewPAIModel(	const G4ParticleDefinition* part, 
                              			const G4String& modname,
                              			const G4String& procname )
	{
  	G4String partname = part->GetParticleName();
  	if(modname == "pai") 
  		{
    		G4PAIModel* pai = new G4PAIModel(part,"PAIModel");
    		fConfig->SetExtraEmModel(partname,procname,pai,"GasDetector", 0.0,100.*TeV, pai);
    	} 
    else if(modname == "pai_photon") 
    	{
    		G4PAIPhotModel* pai = new G4PAIPhotModel(part,"PAIPhotModel");
    		fConfig->SetExtraEmModel(partname,procname,pai,"GasDetector", 0.0,100.*TeV,pai);
  		}
	}

