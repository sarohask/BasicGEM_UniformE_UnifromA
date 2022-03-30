
#include "TrGEMEventAction.hh"
#include "TrGEMAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

TrGEMEventAction::TrGEMEventAction()
	: G4UserEventAction(), DetectorCollID(-1)
	{}

TrGEMEventAction::~TrGEMEventAction()
	{}

void TrGEMEventAction::BeginOfEventAction(const G4Event*)
	{
  	if (DetectorCollID == -1)
  		{
    		DetectorCollID = G4SDManager::GetSDMpointer()->GetCollectionID("TrGEMHitsCollection");
  		}

 	 	TrGEMAnalysis::GetInstance()->PrepareNewEvent() ;
	}

void TrGEMEventAction::EndOfEventAction(const G4Event* event)
	{
		G4int i;
  	G4int numOfDetectHits = 0;
  	
  	//TrGEMAnalysis::GetInstance()->GetNumInteraction();
  	
		G4int event_id = event->GetEventID();
  	if (event_id % 1000000 == 0)
			{
     		printf ("Processed %i events \n", event_id);
    	}
 	
  	TrGEMHitsCollection *trGEMHits = 0;

  	G4HCofThisEvent *HCE = event->GetHCofThisEvent();
		
  	if (HCE)
  		{
    		trGEMHits = (TrGEMHitsCollection *)(HCE->GetHC(DetectorCollID));
    		if (trGEMHits)
    			{
						numOfDetectHits = trGEMHits->GetSize();

      			if (numOfDetectHits > 0)
      				{
      					for (i = 0; i < numOfDetectHits; i++)
      					{
      						TrGEMHit* aHit = (TrGEMHit*)trGEMHits->GetHit(i);
           				TrGEMAnalysis::GetInstance()->SetValuesForInteractingEvents(event_id, aHit);
           				TrGEMAnalysis::GetInstance()->EndOfEvent() ;
           			}		
           		}
      			else	
      				{
      					TrGEMAnalysis::GetInstance()->SetValuesForNonInteractingEvents(event_id);
      					TrGEMAnalysis::GetInstance()->EndOfEvent() ;
        			}
        	}
      }
		
	}
