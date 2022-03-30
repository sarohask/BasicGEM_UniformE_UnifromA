#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT.hh"
// #include "TROOT.h"
#include "G4EmLivermorePhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

#include "TrGEMDetectorConstruction.hh"
#include "TrGEMActionInitialization.hh"
#include "TrGEMAnalysis.hh"
#include "TrGEMPhysicsList.hh"

#include "Randomize.hh"

int main(int argc, char** argv) 
	{
		// Set the Random engine
  		G4Random::setTheEngine(new CLHEP::RanecuEngine());
  		G4Random::setTheSeed(time(NULL)+38999008.);  
  
  		G4RunManager* runManager = new G4RunManager;

  		runManager->SetUserInitialization(new TrGEMDetectorConstruction ) ;
  		G4VUserPhysicsList* physics = new FTFP_BERT_HP();

  		runManager->SetUserInitialization(physics);
  		
  		runManager->SetUserInitialization(new TrGEMActionInitialization(argv[1]));

  		// initialize G4 kernel
  		runManager->Initialize();
  		
   		G4VisManager* visManager = new G4VisExecutive;
   		visManager->Initialize();   
          
  		// get the pointer to the UI manager
  		G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
 		if ( argc == 4 ) 
 			{
 				G4String temp = argv[2];
  				TrGEMAnalysis::GetInstance()->SetFileName(temp);
  				G4String command = "/run/beamOn ";
  				temp = argv[3];
  				//UImanager->ApplyCommand("/run/initialize");
  				//UImanager->ApplyCommand("/cuts/setLowEdge 10 eV");
  				//UImanager->ApplyCommand("/control/verbose 0");
  				//UImanager->ApplyCommand("/run/verbose 0");
  				//UImanager->ApplyCommand("/event/verbose 0");
  				//UImanager->ApplyCommand("/tracking/verbose 2");
  				UImanager->ApplyCommand(command+temp);
			}
		else
			{
		 		// interactive mode : define UI session
				G4UIExecutive* ui = new G4UIExecutive(argc, argv);
	
      			UImanager->ApplyCommand("/control/execute init_vis.mac");     
		
				ui->SessionStart();
				delete ui;
  		
      			delete visManager;  
			}
  		// job termination
  		delete runManager;
  		return 0;
	}
	
	
