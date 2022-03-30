#ifndef TrGEMSteppingAction_h
#define TrGEMSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class TrGEMSteppingAction : public G4UserSteppingAction
	{
   	public:
     	TrGEMSteppingAction();
      ~TrGEMSteppingAction();
			void SetSteppingManagerPointer(G4SteppingManager*);
			void UserSteppingAction(const G4Step*);
			
		private:
			G4SteppingManager* stepManager;
	};

#endif
