
#ifndef TrGEMPhysicsList_h
#define TrGEMPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4EmConfigurator.hh"

class G4VPhysicsConstructor;
class StepMax;
class TrGEMPhysicsListMessenger;

class TrGEMPhysicsList: public G4VModularPhysicsList
	{
		public:
			TrGEMPhysicsList();
  		virtual ~TrGEMPhysicsList();

  		void ConstructParticle();  
  		void SetCuts();
  		void SetCutForGamma(G4double);
  		void SetCutForElectron(G4double);
  		void SetCutForPositron(G4double);
  		void SetCutForProton(G4double);  
  		void AddTrGEMPhysicsList(const G4String& name);
  		void ConstructProcess();
   		void AddStepMax(); 

		private:
			void  AddPAIModel(const G4String&);
  		void  NewPAIModel(const G4ParticleDefinition*, const G4String& modname, const G4String& procname);

  		G4EmConfigurator* fConfig;

			G4double fCutForGamma;
  		G4double fCutForElectron;
  		G4double fCutForPositron;
  		G4double fCutForProton;

  		G4VPhysicsConstructor*  fEmTrGEMPhysicsList;
  		G4VPhysicsConstructor*  fDecayTrGEMPhysicsList;
  		std::vector<G4VPhysicsConstructor*> fHadronPhys;
  		G4String fEmName;
    
  		StepMax* fStepMaxProcess;
	};

#endif

