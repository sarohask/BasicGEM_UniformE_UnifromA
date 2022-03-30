#ifndef TrGEMRunAction_h
#define TrGEMRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class TrGEMRunAction : public G4UserRunAction
	{
  	public:
    	TrGEMRunAction();
    	virtual ~TrGEMRunAction();

  	public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
	};

#endif





