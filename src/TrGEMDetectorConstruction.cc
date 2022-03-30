
#include "TrGEMDetectorConstruction.hh"
#include "TrGEMSensitiveDetector.hh"

#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Para.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "G4EqMagElectricField.hh"
#include "G4SystemOfUnits.hh"

#include "G4UniformElectricField.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"

TrGEMDetectorConstruction::TrGEMDetectorConstruction() 
	: fFR4Mat(0), fGasMat(0), fEmptyMat(0), fAirMat(0), fCuMat(0),fKaptonMat(0), fCoverMat(0), fCoolMat(0), fChipMat(0), fPullMat(0)    
	{
  	driftThinBase  = 282.174*mm ;					 						
  	driftLargeBase = 509.984*mm ;					
  	driftHeight    = 1283.0*mm ;																				
	}

TrGEMDetectorConstruction::~TrGEMDetectorConstruction() 
	{}

void TrGEMDetectorConstruction::DefineMaterials() 
	{
  	G4NistManager* manager = G4NistManager::Instance() ;
  	// define Elements
  	G4Element* elC  = manager->FindOrBuildElement(6);
  	G4Element* elF  = manager->FindOrBuildElement(9);
  	G4Element* elSi = manager->FindOrBuildElement(14);
  	G4Element* elO  = manager->FindOrBuildElement(8);
  	G4Element* elH  = manager->FindOrBuildElement(1);

  	// define Materials
 	 	G4Material *Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si") ;
  	G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu") ;
  	G4Material *Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al") ;
 	 	G4Material *Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
 	 	G4Material *Water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
 	 	G4Material* Air  = manager->FindOrBuildMaterial("G4_AIR");
 	 	G4Material* Steel  = manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  	
  	G4int nElement(0), nAtoms(0) ;
  	G4double density(0.), temperature(0.), pressure(0.), fractionMass(0.), mixtureDensity(0.0)  ;
  	G4String name, symbol ;
  	
  	//Quartz
  	G4Material* SiO2 =  new G4Material("quartz",density= 2.200*g/cm3, nElement=2);
  	SiO2->AddElement(elSi, nAtoms=1);
  	SiO2->AddElement(elO , nAtoms=2);

  	//from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  	//Epoxy (for FR4 )
  	G4Material* Epoxy = new G4Material("Epoxy" , density = 1.2*g/cm3, nElement=2);
  	Epoxy->AddElement(elH, nAtoms=2);
  	Epoxy->AddElement(elC, nAtoms=2);
      
  	//FR4 (Glass + Epoxy)
  	G4Material* FR4 = new G4Material("FR4"  , density = 1.86*g/cm3, nElement=2);
  	FR4->AddMaterial(Epoxy, fractionMass=0.472);
  	FR4->AddMaterial(SiO2, fractionMass=0.528);
  	
  	// gases at STP conditions 
  	G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  	G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  	G4Material* Empty = manager->FindOrBuildMaterial("G4_Galactic");

  	// CF4 must be defined by hand
  	G4Material* CF4 = new G4Material(name="CF4", density=0.003884*g/cm3, nElement=2, kStateGas, temperature = 273.15*kelvin, pressure=1.0*atmosphere);
 	 	CF4->AddElement(elC, 1) ;
  	CF4->AddElement(elF, 4) ; 

  	// Ar:CO2 (70:30) @ STP conditions
  	mixtureDensity = (Argon->GetDensity() * 70/100.0 + CarbonDioxide->GetDensity() * 30/100.0) ;
  	G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  	ArCO2->AddMaterial(Argon, 0.7) ;
  	ArCO2->AddMaterial(CarbonDioxide, 0.3) ;

  	// Ar:CO2:CF4 (45:15:40) @ STP conditions
  	mixtureDensity = (Argon->GetDensity() * 45/100.0 + CarbonDioxide->GetDensity() * 15/100.0 + CF4->GetDensity() * 40/100.0) ;
  	G4Material *ArCO2CF4 = new G4Material("Ar/CO2/CF4",mixtureDensity,3) ;
  	ArCO2CF4->AddMaterial(Argon, 0.45) ;
  	ArCO2CF4->AddMaterial(CarbonDioxide,0.15) ;
  	ArCO2CF4->AddMaterial(CF4,0.40) ;

  	// Choice of the gas
  	fEmptyMat = Empty ;
  	fCuMat = Cu;
  	fKaptonMat = Kapton;
  	fGasMat = ArCO2;
  	fAirMat = Air ;
  	fFR4Mat = FR4;
  	fCoverMat = Al ;
  	fChipMat = Si ;
  	fCoolMat = Water ;
  	fPullMat = Steel ;
	}

G4VPhysicalVolume* TrGEMDetectorConstruction::Construct() 
	{
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	
  	// Define all materials and set global variables
  	DefineMaterials() ;
  	
  	// Visual attributes
  	G4VisAttributes *driftAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ; 										//Green
  	driftAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *gasAttributes = new G4VisAttributes(G4Colour(1., 0.1, 0.9)) ;										//Magenta
  	gasAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *spacerAttributes = new G4VisAttributes(G4Colour(1., 1., 0.)) ;										//Magenta
  	spacerAttributes->SetForceSolid(true) ;
  	G4VisAttributes *gemAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1)) ;										//Light Brown
  	gemAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *kaptonAttributes = new G4VisAttributes(G4Colour(0.7, 0.5, 0.)) ;										//Dark Brown
  	kaptonAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *readoutAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ; 										//Green
  	readoutAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *vFatAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ;
   	vFatAttributes->SetForceSolid(true) ;
   	G4VisAttributes *zigAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)) ;
   	zigAttributes->SetForceSolid(true) ;
   	G4VisAttributes *pullOutAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
   	pullOutAttributes->SetForceSolid(true) ;
   	G4VisAttributes *frameAttributes = new G4VisAttributes(G4Colour(0.3, 0.8, 0.8)) ;
   	frameAttributes->SetForceSolid(true) ;
   	G4VisAttributes *copperAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1));
  	copperAttributes->SetForceWireframe(true);
  	G4VisAttributes *emptyAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0)) ;
   	emptyAttributes->SetForceWireframe(true) ;
   	G4VisAttributes *fakeAttributes = new G4VisAttributes(G4Colour(0., 0., 1.)) ;
  	fakeAttributes->SetForceSolid(true) ;
  	G4VisAttributes *coolingPadAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1)) ;
   	coolingPadAttributes->SetForceSolid(true) ;
  	G4VisAttributes *coolingAttributes = new G4VisAttributes(G4Colour(0.9, 0.6, 0.)) ;
   	coolingAttributes->SetForceWireframe(true) ;
   	G4VisAttributes *coverAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)) ;
  	coverAttributes->SetForceSolid(true) ;
  	G4VisAttributes *worldAttributes = new G4VisAttributes ;
  	worldAttributes->SetVisibility(false) ;
	
  	// rotation Matrix for layers
  	G4RotationMatrix* rotationPlacement = new G4RotationMatrix() ;
  	rotationPlacement->rotateY(M_PI / 2.0) ;
  	rotationPlacement->rotateX(M_PI / 2.0) ;
		
		//roatation matrix for fake lateral / external frame
  	G4RotationMatrix* rotationPlacement1 = new G4RotationMatrix() ;
   	rotationPlacement1->rotateY(M_PI / 2.0) ;
   	rotationPlacement1->rotateX(M_PI / 2.0) ;
   	
   	// Rotation Matrix for bend pipes
   	G4RotationMatrix* rotationPlacement2 = new G4RotationMatrix() ;
  	rotationPlacement2->rotateX(M_PI ) ;
   
   //rotaion matrix for staright pipes
   	G4RotationMatrix* rotationPlacement3 = new G4RotationMatrix() ;
   	rotationPlacement3->rotateY(0.) ;
   	rotationPlacement3->rotateX(M_PI/2.0 ) ;
   	rotationPlacement3->rotateZ(0.) ;
   	
   	//rotation matrix for vfat
   	G4RotationMatrix* rotationPlacement4 = new G4RotationMatrix() ;  	
   	rotationPlacement4->rotateZ(M_PI/2.0 ) ;
   	rotationPlacement4->rotateY(M_PI/2.0) ;
   	
   	//for frame lateral side
   	G4RotationMatrix* rotationPlacement5 = new G4RotationMatrix() ;  	
   	rotationPlacement5->rotateX(M_PI - atan(1283.0/113.905) ) ;
   	
   	//for frame lateral side
   	G4RotationMatrix* rotationPlacement6 = new G4RotationMatrix() ;  	
   	rotationPlacement6->rotateX(-M_PI + atan(1283.0/113.905) ) ;
   	
   	//for pull-outs lateral side
   	G4RotationMatrix* rotationPlacement7 = new G4RotationMatrix() ;  	
   	rotationPlacement7->rotateZ(M_PI - atan(1283.0/113.905) ) ;
   	
   	//for pull-outs lateral side
   	G4RotationMatrix* rotationPlacement8 = new G4RotationMatrix() ;  	
   	rotationPlacement8->rotateZ(-M_PI + atan(1283.0/113.905) ) ;
		
		G4RotationMatrix* rotationPlacement9 = new G4RotationMatrix() ;  	
   	rotationPlacement9->rotateZ(M_PI )  ;
   	
   	G4RotationMatrix* rotationPlacement10 = new G4RotationMatrix() ;  	
   	rotationPlacement10->rotateZ(M_PI+0.959931 )  ;
   	
   	G4RotationMatrix* rotationPlacement11 = new G4RotationMatrix() ;  	
   	rotationPlacement11->rotateZ(0.959931 )  ;
   	
   	G4RotationMatrix* rotationPlacement12 = new G4RotationMatrix() ;  
   	rotationPlacement12->rotateX(M_PI/2. )  ;	
   	rotationPlacement12->rotateY(M_PI/2. -0.60434096913)  ;
   	
   	G4RotationMatrix* rotationPlacement13 = new G4RotationMatrix() ;  
   	rotationPlacement13->rotateX(M_PI/2. )  ;	
   	rotationPlacement13->rotateY(M_PI/2. + 0.60434096913)  ;
   	
		//world
  	G4double worldSizeX = 1.0*m;
  	G4double worldSizeY = 1.0*m;
  	G4double worldSizeZ = 1.0*m;
  	
  	//world definition and placement
  	worldBox = new G4Box("WorldBox", worldSizeX, worldSizeY, worldSizeZ) ;
  	worldLog = new G4LogicalVolume( worldBox, fAirMat, "WorldLog") ;
  	worldLog->SetVisAttributes(worldAttributes) ;
		G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), worldLog, "WorldSpace", 0, false, 0) ;
   
  	std::string layerNames[19]= { "DriftCopper1","DriftBoard","DriftCopper2",     			 	//Drift Board
    															"GasGap1",                                      				//Drift Gap
    															"Gem1Copper1","Gem1","Gem1Copper2",             			 	//GEM1
    															"GasGap2",                                      				//Transfer I Gap
    															"Gem2Copper1","Gem2","Gem2Copper2",             			 	//GEM2
    															"GasGap3",                                      				//Transfer II Gap
    															"Gem3Copper1","Gem3","Gem3Copper2",             			 	//GEM3
    															"GasGap4",                                      				//Induction Gap
    															"ReadCopper1","ReadoutBoard","ReadCopper2",     			 	//Readout Board														 	
  															};
             
  	std::string layerNamesLog[19];
  	for(size_t i=1; i<19; i++) 
  		{ 
    		layerNamesLog[i]=layerNames[i]+"Log";
			}
		
  	G4Material* layerMat[19]= { fCuMat,fFR4Mat,fCuMat,      																//Drift Board
    														fGasMat,                    																//Drift Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM1
    														fGasMat,                    																//Transfer I Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM2
    														fGasMat,                    																//Transfer II Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM3
    														fGasMat,                    																//Induction Gap
    													  fCuMat,fFR4Mat,fCuMat      																//Readout Board
  														};

  	G4double thickness[19] = { 35.*um,3.2*mm,35.*um,      																	//Drift Board
    													 3.*mm,                     																	//Drift Gap
    													 5.*um,50*um,5.*um,         																	//gem1
    													 1.*mm,                     																	//Transfer I Gap
    													 5.*um,50*um,5.*um,         																	//gem2
    													 2.*mm,                     																	//Transfer II Gap
    													 5.*um,50.*um,5.*um,        																	//gem3
    													 1.*mm,                     																	//Induction Gap
    													 35.*um,3.2*mm,35.*um      																	//Readout Board
  													 };
  													
  	G4double layerLargeBase[19] = { driftLargeBase,driftLargeBase,driftLargeBase,      			//Drift Board
    																driftLargeBase,                     											//Drift Gap
    																driftLargeBase,driftLargeBase,driftLargeBase,         		//gem1
    																driftLargeBase,                     											//Transfer I Gap
    																driftLargeBase,driftLargeBase,driftLargeBase,         		//gem2
    																driftLargeBase,                    												//Transfer II Gap
    																driftLargeBase,driftLargeBase,driftLargeBase,        		//gem3
    																driftLargeBase,                     											//Induction Gap
    																driftLargeBase,driftLargeBase,driftLargeBase     		//Readout Board
  																};
  	
  	G4double layerThinBase[19] = { driftThinBase,driftThinBase,driftThinBase,      					//Drift Board
    															 driftThinBase,                     												//Drift Gap
    															 driftThinBase,driftThinBase,driftThinBase,         				//gem1
    															 driftThinBase,                     												//Transfer I Gap
    															 driftThinBase,driftThinBase,driftThinBase,         				//gem2
    															 driftThinBase,                     												//Transfer II Gap
    															 driftThinBase,driftThinBase,driftThinBase,        				//gem3
    															 driftThinBase,                     												//Induction Gap
    															 driftThinBase,driftThinBase,driftThinBase     				//Readout Board
  															 };
  	
  	G4double layerHeight[19] = { driftHeight,driftHeight,driftHeight,      									//Drift Board
    														 driftHeight,                     														//Drift Gap
    														 driftHeight,driftHeight,driftHeight,         								//gem1
    														 driftHeight,                     														//Transfer I Gap
    														 driftHeight,driftHeight,driftHeight,         								//gem2
    														 driftHeight,                     														//Transfer II Gap
    														 driftHeight,driftHeight,driftHeight,        								//gem3
    														 driftHeight,                   															//Induction Gap
    														 driftHeight,driftHeight,driftHeight     								//Readout Board
  														};
  	
  	G4VisAttributes *visAttributes[19] = { copperAttributes, driftAttributes, copperAttributes,					//Drift Board
  																		 		 gasAttributes,																								//Drift Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem1
  																		 		 gasAttributes,																								//Transfer I Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem2
  																		 		 gasAttributes,																								//Transfer II Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem3
  																		 		 gasAttributes,																								//Induction Gap
  																		 		 copperAttributes, readoutAttributes, copperAttributes,				//Readout Board
  																			};
  																
  																			
  	// Chamber values
  	double chamberZPos[numChamber] = { -15.57, 1.85 };
  	
  	// Fake values
  	double fakeThickness = 0.1*mm ;
  	double fakeZPos[numChamber][numFake] = {{-15.62*mm, -1.8*mm}, {1.8*mm, 15.62*mm}} ;
  	
  	// For loop --- 2 Chambers														 	
  	for(G4int chamber = 0; chamber < buildChamber ; chamber++)
  		{				
				G4Trd* strato[numChamber];
  			G4LogicalVolume* logicStrato[numChamber];
  
  			for(G4int lyr=0; lyr<19; lyr++)
  				{
    				strato[chamber] = Trapezoid(layerNames[lyr], thickness[lyr],layerLargeBase[lyr],layerThinBase[lyr], layerHeight[lyr]) ;
    				logicStrato[chamber] = new G4LogicalVolume (strato[chamber], layerMat[lyr],layerNamesLog[lyr]);   
    				logicStrato[chamber]->SetVisAttributes(visAttributes[lyr]) ;
    				trdCollection[chamber].push_back(strato[chamber]) ;
    				trdLogCollection[chamber].push_back(logicStrato[chamber]) ;
  				} 
  		
				PlaceGeometry(chamber, rotationPlacement,G4ThreeVector(0.,0.,chamberZPos[chamber]),worldLog) ;
				
				for(G4int oFake = 0; oFake < 2; oFake++)
  				{				
						fake[chamber][oFake] = new G4Trd("Fake"+to_string(oFake)+to_string(chamber),  fakeThickness/2., fakeThickness/2., driftThinBase/2., driftLargeBase/2., driftHeight/2.);
  					fakeLog[chamber][oFake] = new G4LogicalVolume(fake[chamber][oFake], fAirMat, "FakeLog");
						fakeLog[chamber][oFake] ->SetVisAttributes(new G4VisAttributes(*fakeAttributes)) ;
			
						new G4PVPlacement ( rotationPlacement,
									  						G4ThreeVector(0., 0., fakeZPos[chamber][oFake]),
									  						fakeLog[chamber][oFake] ,
									  						"Fake"+to_string(oFake)+"Phy"+to_string(chamber),
									  						worldLog,
									 		 					false,
									  						0);	
					}
			}
  	return worldPhys ;
	}

G4Trd* TrGEMDetectorConstruction::Trapezoid(G4String name, G4double width, G4double largeBase, G4double thinBase, G4double height) 
	{
  	G4Trd* shape = new G4Trd( name,
                           	 	width/2, width/2,
                           		thinBase/2,
                           		largeBase/2,
                           		height/2) ;
  	return shape ;
	}

void TrGEMDetectorConstruction::PlaceGeometry(G4int chamberNum, G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume* pMotherLogical) 
	{
  	G4double XTranslation = 0. ;
   	
   	for(size_t i=0 ; i<trdCollection[chamberNum].size() ; i++) 
   		{
      	// i counts as the copyNo
      	G4String layerName = trdCollection[chamberNum].at(i)->GetName() ;
      	XTranslation += trdCollection[chamberNum].at(i)->GetXHalfLength1() ;
      	G4ThreeVector position = tlate + G4ThreeVector(XTranslation,0,0).transform(G4RotationMatrix(*pRot).inverse()) ;
      	G4cout <<"For Chamber "<<chamberNum<< ".......Volume (" << i << ") " << layerName << " the position is " << G4BestUnit(XTranslation,"Length") << G4endl ;
  
				new G4PVPlacement(pRot, position, trdLogCollection[chamberNum].at(i), layerName+"Phy"+to_string(chamberNum), pMotherLogical, false, i) ;

   			XTranslation += trdCollection[chamberNum].at(i)->GetXHalfLength1() ;
			}
	}

void TrGEMDetectorConstruction::ConstructSDandField()
	{
		G4SDManager* sdman = G4SDManager::GetSDMpointer() ;
		
		TrGEMSensitiveDetector* sensitive = new TrGEMSensitiveDetector("/GasGap") ;
  	sdman->AddNewDetector(sensitive) ;
  	
  	for(G4int chamber = 0; chamber < buildChamber ; chamber++)
  		{
  			fakeLog[chamber][0]->SetSensitiveDetector(sensitive) ;
  			fakeLog[chamber][1]->SetSensitiveDetector(sensitive) ;
				trdLogCollection[chamber][3]->SetSensitiveDetector(sensitive) ;									//GasGap1Log	
  			trdLogCollection[chamber][7]->SetSensitiveDetector(sensitive) ;									//GasGap2Log	
  			trdLogCollection[chamber][11]->SetSensitiveDetector(sensitive) ;								//GasGap3Log
  			trdLogCollection[chamber][15]->SetSensitiveDetector(sensitive) ;								//GasGap4Log		
 		}
	}
