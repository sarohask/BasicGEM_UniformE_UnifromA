#ifndef TGEMTrackingAction_h
#define TrGEMTrackingAction_h 1

#include "globals.hh"
#include "G4UserTrackingAction.hh"

using namespace std;

class TrGEMTrackingAction : public G4UserTrackingAction
	{
   	public:
      TrGEMTrackingAction() ;
      virtual ~TrGEMTrackingAction() ;

      virtual void PreUserTrackingAction(const G4Track*) ;
      virtual void PostUserTrackingAction(const G4Track*) ;
	};

#endif
