
#include "GEMDPG/TrGEMG4/interface/TrGEMDetectorConstruction.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMSensitiveDetector.hh"

#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
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
	: fFR4Mat(0), fGasMat(0), fEmptyMat(0), fAirMat(0), fCuMat(0),fKaptonMat(0), tripleGemThinBase(0), tripleGemLargeBase(0), tripleGemHeight(0)    
	{
  	tripleGemThinBase  = 279.0*mm ;
  	tripleGemLargeBase = 510.0*mm ;
  	tripleGemHeight    = 1283.0*mm ;
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
  	G4Material *H  = G4NistManager::Instance()->FindOrBuildMaterial("G4_H") ;
  	G4Material *C  = G4NistManager::Instance()->FindOrBuildMaterial("G4_C") ;
  	G4Material *O  = G4NistManager::Instance()->FindOrBuildMaterial("G4_O") ;
 	 	G4Material *Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si") ;
  	G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu") ;
 	 	G4Material *Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
 	 	G4Material* Air  = manager->FindOrBuildMaterial("G4_AIR");
  	
  	G4int numel(0), natoms(0) ;
  	G4double density(0.), temperature(0.), pressure(0.), fractionMass(0.), mixtureDensity(0.0)  ;
  	G4String name, symbol ;
  	
  	//Quartz
  	G4Material* SiO2 =  new G4Material("quartz",density= 2.200*g/cm3, numel=2);
  	SiO2->AddElement(elSi, natoms=1);
  	SiO2->AddElement(elO , natoms=2);

  	//from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  	//Epoxy (for FR4 )
  	G4Material* Epoxy = new G4Material("Epoxy" , density = 1.2*g/cm3, numel=2);
  	Epoxy->AddElement(elH, natoms=2);
  	Epoxy->AddElement(elC, natoms=2);
      
  	//FR4 (Glass + Epoxy)
  	G4Material* FR4 = new G4Material("FR4"  , density = 1.86*g/cm3, numel=2);
  	FR4->AddMaterial(Epoxy, fractionMass=0.472);
  	FR4->AddMaterial(SiO2, fractionMass=0.528);
  	
  	// gases at STP conditions 
  	G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  	G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  	G4Material* empty = manager->FindOrBuildMaterial("G4_Galactic");
  	fEmptyMat = empty ;

  	// CF4 must be defined by hand
  	G4Material* CF4 = new G4Material(name="CF4", density=0.003884*g/cm3, numel=2, kStateGas, temperature = 273.15*kelvin, pressure=1.0*atmosphere);
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
  	fCuMat = Cu;
  	fKaptonMat = Kapton;
  	fGasMat = ArCO2;
  	fAirMat = Air ;
  	fFR4Mat = FR4;
	}

G4VPhysicalVolume* TrGEMDetectorConstruction::Construct() 
	{
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	
  	// Define all materials and set global variables
  	DefineMaterials() ;

  	// SD Manager 
  	G4SDManager* sdman = G4SDManager::GetSDMpointer() ;

  	G4double worldSizeX = 1.0*m;
  	G4double worldSizeY = 1.0*m;
  	G4double worldSizeZ = 1.0*m;

  	// Rotation Matrix for layers
  	G4RotationMatrix* rotationPlacement = new G4RotationMatrix() ;
  	rotationPlacement->rotateY(M_PI / 2.0) ;
  	rotationPlacement->rotateX(M_PI / 2.0) ;

  	// Visual attributes
  	G4VisAttributes *cathodeAttributes = new G4VisAttributes(G4Color::Grey()) ;
  	cathodeAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *g10Attributes = new G4VisAttributes(G4Color::White()) ;
  	g10Attributes->SetForceWireframe(true) ;
  	G4VisAttributes *gasAttributes = new G4VisAttributes(G4Color::Red()) ;
  	gasAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *gemAttributes = new G4VisAttributes(G4Color::Green()) ;
  	gemAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *insAttributes = new G4VisAttributes(G4Color::Blue()) ;
  	insAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *worldAttributes = new G4VisAttributes ;
  	worldAttributes->SetVisibility(true) ;

  	// World definition and placement
  	G4Box* worldBox 							= new G4Box("WorldBox", worldSizeX, worldSizeY, worldSizeZ) ;
  	G4LogicalVolume* worldLog 		= new G4LogicalVolume( worldBox, fEmptyMat, "WorldLog") ;
  		worldLog->SetVisAttributes(worldAttributes) ;
		G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), worldLog, "WorldSpace", 0, false, 0) ;
   
  	std::string NomeStrati[21]= { "FakeBottom",                                   //Fake
    															"DriftCopper1","DriftBoard","DriftCopper2",     //Drift Board
    															"GasGap1",                                      //Drift Gap
    															"Gem1Copper1","Gem1","Gem1Copper2",             //gem1
    															"GasGap2",                                      //Transfer I Gap
    															"Gem2Copper1","Gem2","Gem2Copper2",             //gem2
    															"GasGap3",                                      //Transfer II Gap
    															"Gem3Copper1","Gem3","Gem3Copper2",             //gem3
    															"GasGap4",                                      //Induction Gap
    															"ReadCopper1","ReadoutBoard","ReadCopper2",     //Readout Board
    															"FakeTop"                                       //Fake
  															};
             
  	std::string NomeStratiLog[21];
  	for(size_t i=1; i<21; i++) 
  		{ 
    		NomeStratiLog[i]=NomeStrati[i]+"Log";
			}
		
  	G4Material* MatStrati[21]= { fAirMat,                    //Fake
    														 fCuMat,fFR4Mat,fCuMat,      //Drift Board
    														 fGasMat,                    //Drift Gap
    														 fCuMat,fKaptonMat,fCuMat,   //gem1
    														 fGasMat,                    //Transfer I Gap
    														 fCuMat,fKaptonMat,fCuMat,   //gem2
    														 fGasMat,                    //Transfer II Gap
    														 fCuMat,fKaptonMat,fCuMat,   //gem3
    														 fGasMat,                    //Induction Gap
    													   fCuMat,fFR4Mat,fCuMat,      //Readout Board
    												     fAirMat                     //Fake
  														};

  	G4double spessoreStrati[21] = { 0.1*mm,                    //Fake
    																35.*um,3.2*mm,35.*um,      //Drift Board
    																3.*mm,                     //Drift Gap
    																5.*um,50*um,5.*um,         //gem1
    																1.*mm,                     //Transfer I Gap
    																5.*um,50*um,5.*um,         //gem2
    																2.*mm,                     //Transfer II Gap
    																5.*um,50.*um,5.*um,        //gem3
    																1.*mm,                     //Induction Gap
    																35.*um,3.2*mm,35.*um,      //Readout Board
    																0.1*mm                     //Fake
  																};

		TrGEMSensitiveDetector* sensitive = new TrGEMSensitiveDetector("/GasGap") ;
  	sdman->AddNewDetector(sensitive) ;
  
  	G4Trd* strato;
  	G4LogicalVolume* logicStrato;
  	
  	for(G4int lyr=0; lyr<21; lyr++)
  		{
    		strato = Trapezoid(NomeStrati[lyr], spessoreStrati[lyr]);
    		logicStrato = new G4LogicalVolume (strato, MatStrati[lyr],NomeStratiLog[lyr]);   
    		logicStrato->SetVisAttributes(new G4VisAttributes(*gemAttributes)) ;
    		trdCollection.push_back(strato) ;
    		trdLogCollection.push_back(logicStrato) ;
  		}
  	trdLogCollection[0]->SetSensitiveDetector(sensitive) ;	
  	//trdLogCollection[1]->SetSensitiveDetector(sensitive) ;	
  	//trdLogCollection[2]->SetSensitiveDetector(sensitive) ;
  	//trdLogCollection[3]->SetSensitiveDetector(sensitive) ;
  	trdLogCollection[4]->SetSensitiveDetector(sensitive) ;								//GasGap1Log	
  	//trdLogCollection[5]->SetSensitiveDetector(sensitive) ;	
  	//trdLogCollection[6]->SetSensitiveDetector(sensitive) ;	
  	//trdLogCollection[7]->SetSensitiveDetector(sensitive) ;		
  	trdLogCollection[8]->SetSensitiveDetector(sensitive) ;								//GasGap2Log   
  	//trdLogCollection[9]->SetSensitiveDetector(sensitive) ;		   
  	//trdLogCollection[10]->SetSensitiveDetector(sensitive) ;		   
  	//trdLogCollection[11]->SetSensitiveDetector(sensitive) ;	
  	trdLogCollection[12]->SetSensitiveDetector(sensitive) ;								//GasGap3Log
  	//trdLogCollection[13]->SetSensitiveDetector(sensitive) ;
  	//trdLogCollection[14]->SetSensitiveDetector(sensitive) ;
  	//trdLogCollection[15]->SetSensitiveDetector(sensitive) ;
  	trdLogCollection[16]->SetSensitiveDetector(sensitive) ;								//GasGap4Log
  	//trdLogCollection[17]->SetSensitiveDetector(sensitive) ;
  	//trdLogCollection[18]->SetSensitiveDetector(sensitive) ;
  	//trdLogCollection[19]->SetSensitiveDetector(sensitive) ;
  	trdLogCollection[20]->SetSensitiveDetector(sensitive) ; 
  
		PlaceGeometry(rotationPlacement,G4ThreeVector(0.,0.,0.),worldLog) ;
  	return worldPhys ;
	}

G4Trd* TrGEMDetectorConstruction::Trapezoid(G4String name, G4double width) 
	{
  	G4Trd* shape = new G4Trd( name,
                           	 	width/2, width/2,
                           		tripleGemThinBase/2,
                           		tripleGemLargeBase/2,
                           		tripleGemHeight/2) ;
  	return shape ;
	}

void TrGEMDetectorConstruction::PlaceGeometry(G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume* pMotherLogical) 
	{
  	G4double XTranslation = 0.0 ;
  	for(size_t i=0 ; i<trdCollection.size() ; i++) 
  		{
    		// i counts as the copyNo
    		XTranslation += trdCollection.at(i)->GetXHalfLength1() ;
    		G4ThreeVector position = tlate + G4ThreeVector(XTranslation,0,0).transform(G4RotationMatrix(*pRot).inverse()) ;
    		G4cout << "Volume (" << i << ") " << trdCollection.at(i)->GetName() << " the position is " << G4BestUnit(XTranslation,"Length") << G4endl ;
    		new G4PVPlacement(pRot,
                      		position,
                      		trdLogCollection.at(i),
                      		trdCollection.at(i)->GetName(),
                      		pMotherLogical,
                      		false,
                      		i) ;
    		XTranslation += trdCollection.at(i)->GetXHalfLength1() ;
  		}
	}

