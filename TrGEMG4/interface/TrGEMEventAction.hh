#ifndef TrGEMEventAction_h
#define TrGEMEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh" 

class G4Event;

class TrGEMEventAction : public G4UserEventAction
	{
   	public:
      TrGEMEventAction();
      ~TrGEMEventAction();
      void BeginOfEventAction(const G4Event* anEvent);
      void EndOfEventAction(const G4Event* anEvent);
      
  	private:
   		G4int DetectorCollID ;
	};

#endif


