///////////////////////////////////////////////////////////////////
//
// Apr/2015  E. Nacher -> DetectorConstruction.cc
//
///////////////////////////////////////////////////////////////////

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

using namespace CLHEP;
namespace B1
{

double cross_product(G4ThreeVector* a, G4ThreeVector* b){
	return a->getX()*b->getX()+a->getY()*b->getY()+a->getZ()*b->getZ();
}

double magnitude(G4ThreeVector* a){
	return sqrt(a->getX()*a->getX()+a->getY()*a->getY()+a->getZ()*a->getZ());
}

G4ThreeVector* vector_product(G4ThreeVector* a, G4ThreeVector* b){
	double X, Y, Z;
	X = a->getY()*b->getZ() - a->getZ()*b->getY();
	Y = a->getZ()*b->getX() - a->getX()*b->getZ();
	Z = a->getX()*b->getY() - a->getY()*b->getX();
	auto vec = new G4ThreeVector(X,Y,Z);
	return vec;
}

double angle_between(G4ThreeVector* a, G4ThreeVector* b){
	double up = cross_product(a,b);
	G4cout << "UP: " << up << G4endl;
	double down = magnitude(a)*magnitude(b);
	return acos(up/down);//*180/M_PI;
}

/*DetectorConstruction::DetectorConstruction()
{ }

DetectorConstruction::~DetectorConstruction()
{ }
*/

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	G4NistManager* nist = G4NistManager::Instance();
	G4RotationMatrix* rotm_90_deg = new G4RotationMatrix(0, CLHEP::pi/2, 0);
	G4bool checkOverlaps = true;
	
	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------
	
	G4String name, symbol;             //a=mass of a mole;
	G4double a, z, density;            //z=mean number of protons;  
	
	G4int ncomponents, natoms;         
	
	G4double pressure    = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	density     = 1.e-25*g/cm3;
	
	G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, 
										density,kStateGas,temperature,pressure);
	
	
	//
	// define simple Elements
	//
	
	// O
	a = 15.999*g/mole;
	G4Element* elO  = new G4Element(name="Oxygen" ,symbol="O" , z= 8., a);
		
	// Mg
	a = 24.305*g/mole;
	G4Element* elMg  = new G4Element(name="Magnesium"  ,symbol="Mg" , z= 12., a);
	
	// Si
	a = 28.0855*g/mole;
	G4Element* elSi  = new G4Element(name="Silicon"  ,symbol="Si" , z= 14., a);
	
	// Cl
	a = 35.453*g/mole;
	G4Element* elCl  = new G4Element(name="Chlorine"  ,symbol="Cl" , z= 17., a);
	
	// K
	a = 39.0983*g/mole;
	G4Element* elK  = new G4Element(name="Potassium"  ,symbol="K" , z= 19., a);
	
	// Br
	a = 79.904*g/mole;
	G4Element* elBr  = new G4Element(name="Bromine"  ,symbol="Br" , z= 35., a);
	
	// Sb
	a = 121.76*g/mole;
	G4Element* elSb = new G4Element(name="Antimony",symbol="Sb" , z= 51., a);
	
	// Cs
	a = 132.905*g/mole;
	G4Element* elCs = new G4Element(name="Cesium",symbol="Cs" , z= 55., a);
	
	// La
	a = 138.905*g/mole;
	G4Element* elLa = new G4Element(name="Lanthanum",symbol="La" , z= 57., a);
	
	G4Element* elC = nist->FindOrBuildElement("C");
	G4Element* elH = nist->FindOrBuildElement("H");
	G4Element* elN = nist->FindOrBuildElement("N");
	
	
	//
	// define simple materials
	//
	
	// Outer sphere material
	G4Material* outer_sphere_mat = nist->FindOrBuildMaterial("G4_Al");
	
	// Cylinder material
	G4Material* cyl_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	
	// Si pmt material
	density = 3.0*g/cm3;
	G4Material* Si_mat = new G4Material(name="Si_pmt_material", density, ncomponents=1);
	Si_mat->AddElement(elSi, natoms=1);
	
	// Out sphere material
	G4Material* out_sphere_mat = nist->FindOrBuildMaterial("G4_Al");
	G4Material* pmt_mat = nist->FindOrBuildMaterial("G4_Glass");

	// Scintillator 1
	G4Material* scint_mat = new G4Material("lin_alkyl_ppo", 0.863*g/cm3, 2); // C6H5C10H21 + PPO
	scint_mat->AddElementByNumberOfAtoms(elC, 15);
	scint_mat->AddElementByNumberOfAtoms(elH, 11);

	// Scintillator 2
	G4Material* lab = new G4Material("lin_alkyl", 0.863*g/cm3, 2);
	lab->AddElementByNumberOfAtoms(elC, 15);
	lab->AddElementByNumberOfAtoms(elH, 11);
	
	// Borosilicate Glass
	G4Material* BSG = new G4Material("Bor_silc_glass", 2.23*g/cm3, 2); 
	BSG->AddElementByNumberOfAtoms(elSi, 1);
	BSG->AddElementByNumberOfAtoms(elO, 2);

	// Glass
	G4Material* PMMA = new G4Material("PMMA", 1.18*g/cm3, 3); // PMMA C5H8O2
	PMMA->AddElementByNumberOfAtoms(elC, 5);
	PMMA->AddElementByNumberOfAtoms(elH, 8);
	PMMA->AddElementByNumberOfAtoms(elO, 2);
	
	//------------------------------------------------------
	// Optical properties
	//------------------------------------------------------

	// Scintillator 1 properties
	
	
	/*
	
	std::vector<G4double> lxe_ABSL  = { 35. * cm, 35. * cm, 35. * cm };
	
	G4double Ephoton_sc[] = {1.0*eV, 7.0*eV};
	G4double refractive_Scint[] = {1.5, 1.5};
	G4double absorption_Scint[] = {5*m, 5*m}; // linalkyl
	//G4double absorption_Scint[] = {10*cm, 10*cm};
	
	
	//std::vector<G4double> lxe_Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };
	
	*/
	
	//G4double Ephoton_sc[] = {1.0*eV, 7.0*eV};
	G4double Ephoton_sc[] = {1.0*eV, 7.0*eV};
	G4double refractive_Scint[] = {1.5, 1.5};
	G4double reflectivity_index[] = {1., 1.};
	
	G4double refractiveIndexPMMA[] =
	{1.49, 1.49};
	
	G4double refractiveIndexBSG[] =
	{1.471, 1.441};
	
	G4double absorptionPMMA[] =
	{10*m, 10*m};
	
	G4double absorptionBSG[] =
	{10*m, 10*m};
	
	G4double absorption_Scint_lab[] = {12*m, 12*m}; // lab
	G4double absorption_Scint[] = {5*m, 5*m}; // linalkyl
	
	G4int nEntries = 12;
	G4double ScintFast[nEntries] = {0,0.03,0.17,0.40,0.55,0.83,1.00,0.84,0.49,0.20,0.07,0.04};
	G4double PhotonEnergy[nEntries] = {2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV, 2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV, 3.26*eV, 3.44*eV};
	
	G4MaterialPropertiesTable* PMMA_mpt = new G4MaterialPropertiesTable();
	PMMA->SetMaterialPropertiesTable(PMMA_mpt);
	PMMA_mpt ->AddProperty("RINDEX", Ephoton_sc, refractiveIndexPMMA, 2);
	PMMA_mpt->AddProperty("ABSLENGTH", Ephoton_sc, absorptionPMMA, 2);
	
	G4MaterialPropertiesTable* BSG_mpt = new G4MaterialPropertiesTable();
	BSG->SetMaterialPropertiesTable(BSG_mpt);
	BSG_mpt ->AddProperty("RINDEX", Ephoton_sc, refractiveIndexBSG, 2);
	BSG_mpt->AddProperty("ABSLENGTH", Ephoton_sc, absorptionBSG, 2);
	
	G4MaterialPropertiesTable* scint_mat_mpt = new G4MaterialPropertiesTable();
	scint_mat_mpt->AddProperty("SCINTILLATIONCOMPONENT1",PhotonEnergy,ScintFast,nEntries);
	scint_mat_mpt->AddConstProperty("SCINTILLATIONYIELD",10000./MeV); 
	scint_mat_mpt->AddConstProperty("RESOLUTIONSCALE",1.);
	scint_mat_mpt->AddConstProperty("SCINTILLATIONRISETIME1", 0.9 * ns);
	scint_mat_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",2.*ns);
	scint_mat_mpt->AddConstProperty("SCINTILLATIONYIELD1",1.);
	scint_mat_mpt->AddConstProperty("SCINTILLATIONYIELD2",0.);
	//scint_mat_mpt->AddProperty("REFLECTIVITY",Ephoton_sc,reflectivity_index,2);

	
  
  
	//scint_mat_mpt -> AddConstProperty("SCINTILLATIONYIELD", 5000./MeV);
	//scint_mat_mpt -> AddConstProperty("ABSLENGTH", 5*m);
	scint_mat_mpt->AddProperty("RINDEX", Ephoton_sc, refractive_Scint, 2);
	scint_mat_mpt->AddProperty("ABSLENGTH", Ephoton_sc, absorption_Scint, 2);
	//scint_mat_mpt -> AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
	//scint_mat->SetMaterialPropertiesTable(scint_mat_mpt);
	
	scint_mat->SetMaterialPropertiesTable(scint_mat_mpt);

	// Scintillator 2 properties
	G4MaterialPropertiesTable* lab_mpt = new G4MaterialPropertiesTable();		 
	lab_mpt -> AddProperty("ABSLENGTH", Ephoton_sc, absorption_Scint_lab, 2);
	lab_mpt ->AddProperty("RINDEX", Ephoton_sc, refractiveIndexPMMA, 2);
	
	lab_mpt->AddProperty("SCINTILLATIONCOMPONENT1",PhotonEnergy,ScintFast,nEntries);
	lab->SetMaterialPropertiesTable(lab_mpt);

	
	//------------------------------------------------------
	// Detector geometry
	//------------------------------------------------------
	
	//     
	// World
	//
	
	G4double WorldSize= 300.*cm;
	
	G4Box* 
    solidWorld = new G4Box("World",		       	                  //its name
						   WorldSize/2,WorldSize/2,WorldSize/2);  //its size
	
	G4LogicalVolume* 
    logicWorld = new G4LogicalVolume(solidWorld,      	//its solid
									 Vacuum,	        //its material
									 "World");		    //its name
	
	G4VPhysicalVolume* 
    physiWorld = new G4PVPlacement(0,			    //no rotation
								   G4ThreeVector(),	//at (0,0,0)
								   "World",		    //its name
								   logicWorld,		//its logical volume
								   NULL,		    //its mother  volume
								   false,	       	//no boolean operation
								   0);			    //copy number
	
	
	//
	// Detector 
	//
	
	//
	// PMT Coordinates
	//
	
	//double f = (1+sqrt(5))/2;
	/*
	double r = 800;
	double pmt_20[20][3] = {
                            {r/sqrt(3), r/sqrt(3), r/sqrt(3)},
                            {-r/sqrt(3), r/sqrt(3), r/sqrt(3)},
                            {-r/sqrt(3), -r/sqrt(3), r/sqrt(3)},
                            {-r/sqrt(3), -r/sqrt(3), -r/sqrt(3)},
                            {r/sqrt(3), -r/sqrt(3), -r/sqrt(3)},
                            {r/sqrt(3), r/sqrt(3), -r/sqrt(3)},
                            {-r/sqrt(3), r/sqrt(3), -r/sqrt(3)},
                            {r/sqrt(3), -r/sqrt(3), r/sqrt(3)},

                            {0, r/(sqrt(3)*f), r*f/sqrt(3)},
                            {0, -r/(sqrt(3)*f), r*f/sqrt(3)},
                            {0, -r/(sqrt(3)*f), -r*f/sqrt(3)},
                            {0, r/(sqrt(3)*f), -r*f/sqrt(3)},

                            {r/(sqrt(3)*f), r*f/sqrt(3), 0},
                            {-r/(sqrt(3)*f), r*f/sqrt(3), 0},
                            {-r/(sqrt(3)*f), -r*f/sqrt(3), 0},
                            {r/(sqrt(3)*f), -r*f/sqrt(3), 0},

                            {r*f/sqrt(3), 0, r/(sqrt(3)*f)},
                            {-r*f/sqrt(3), 0, r/(sqrt(3)*f)},
                            {-r*f/sqrt(3), 0, -r/(sqrt(3)*f)},
                            {r*f/sqrt(3), 0, -r/(sqrt(3)*f)}
                        };
        */
        /*               
        double r_12 = 460 // sqrt(3) factor
        double pmt_12[12][3] = {
			    {f*r,r,0},
			    {f*r,-r,0},
			    {-f*r,-r,0},
			    {-f*r,r,0},
			    
			    {r,0,f*r},
			    {-r,0,f*r},
			    {-r,0,-f*r},
			    {r,0,-f*r},
			    
			    {0,f*r,r},
			    {0,f*r,-r},
			    {0,-f*r,-r},
			    {0,-f*r,r}
			};
	*/
	G4double cylinderHalfHeight = 1860*mm/2;
	G4double cylinderRadius = 1860*mm/2;
	G4double thickness = 2*mm;
	
	G4double sphereRadius = 630.35*mm;
	G4double sphereThickness = 10*mm;
	
	
	//G4double ScintHalfLength =1.5*cm;
	//G4double ScintRadius = 2.*cm;
	
        //G4double ReflectorThickness = 0*mm;
	//G4double ReflectorHalfLength = ScintHalfLength+ReflectorThickness;
	//G4double ReflectorRadius = ScintRadius+ReflectorThickness;
	
	G4double PMTWindowHalfLength = 83.0*mm;
	G4double PMTWindowRadius = 101*mm;
	
	//G4double CathodeHalfLength = 0.005*mm;
	//G4double CathodeRadius =1.9*cm;
	
	G4double StartPhi = 0.*deg;
	G4double DeltaPhi = 360.*deg;
	
	
	// Cylinder
	
	G4Tubs *inner_cylinder = new G4Tubs("inner_cylinder",                    
						0*mm, cylinderRadius-thickness, cylinderHalfHeight-thickness, StartPhi, DeltaPhi); 

	G4LogicalVolume* logic_inner_cylinder = new G4LogicalVolume(inner_cylinder,  
						scint_mat, "logic_inner_cylinder");

	G4VPhysicalVolume* physic_inner_cylinder = new G4PVPlacement(rotm_90_deg,  // yes rotation
								G4ThreeVector(),          // at (0,0,0)
								logic_inner_cylinder,                 // its logical volume
								"inner_cylinder",               // its name
								logicWorld,               // its mother  volume
								false,                    // no boolean operation
								0,                        // copy number
								checkOverlaps);           // overlaps checking

	G4Tubs* outer_cylinder = new G4Tubs("outer_cylinder",                    
						cylinderRadius-thickness, cylinderRadius, cylinderHalfHeight-thickness, StartPhi, DeltaPhi);  

	G4LogicalVolume* logic_outer_cylinder = new G4LogicalVolume(outer_cylinder, 
						cyl_mat, "logic_outer_cylinder");

	G4VPhysicalVolume* physic_outer_cylinder = new G4PVPlacement(rotm_90_deg,  // yes rotation
								G4ThreeVector(),          // at (0,0,0)
								logic_outer_cylinder,                 // its logical volume
								"out_cylinder",               // its name
								logicWorld,               // its mother  volume
								false,                    // no boolean operation
								0,                        // copy number
								checkOverlaps);           // overlaps checking

	// Cylinder lids 
	
	G4Tubs *top_lid = new G4Tubs("top_lid", 0*mm, cylinderRadius, thickness, StartPhi, DeltaPhi);
	
	G4LogicalVolume* logic_top_lid = new G4LogicalVolume(top_lid, 
						cyl_mat, "logic_top_lid"); 

	G4VPhysicalVolume* physic_top_lid = new G4PVPlacement(rotm_90_deg,  // no rotation
						G4ThreeVector(0, 930*mm, 0),     
						logic_top_lid,                 // its logical volume
						"top_lid",               // its name
						logicWorld,               // its mother  volume
						false,                    // no boolean operation
						0,                        // copy number
						checkOverlaps);           // overlaps checking

	G4Tubs *bottom_lid = new G4Tubs("bottom_lid", 0*mm, cylinderRadius, thickness, StartPhi, DeltaPhi);
	
	G4LogicalVolume* logic_bottom_lid = new G4LogicalVolume(bottom_lid, 
						cyl_mat, "logic_bottom_lid"); 

	G4VPhysicalVolume* physic_bottom_lid = new G4PVPlacement(rotm_90_deg,  // no rotation
						G4ThreeVector(0, -930*mm, 0),         
						logic_bottom_lid,                 // its logical volume
						"bottom_lid",               // its name
						logicWorld,               // its mother  volume
						false,                    // no boolean operation
						0,                        // copy number
						checkOverlaps);           // overlaps checking
						
	// Sphere filled with liquid scintillator
	
	// Enveloping part
	
	G4Sphere *outer_sphere = new G4Sphere("outer_sphere", sphereRadius-sphereThickness, sphereRadius, StartPhi, DeltaPhi, StartPhi, DeltaPhi/2); 

	G4LogicalVolume* logic_outer_sphere = new G4LogicalVolume(outer_sphere, 
							PMMA, "outer_sphere"); 

	G4VPhysicalVolume* physic_outer_sphere = new G4PVPlacement(nullptr,  // no rotation
						G4ThreeVector(),                     // at position
						logic_outer_sphere,              // its logical volume
						"outer_sphere",                 // its name
						logic_inner_cylinder,                 // its mother  volume
						true,                    // yes boolean operation
						0,                        // copy number
						checkOverlaps);           // overlaps checking

	// Actual scintillator inside
	
	G4Sphere *inner_sphere = new G4Sphere("inner_sphere", 0.*mm, sphereRadius-sphereThickness, StartPhi, DeltaPhi, StartPhi, DeltaPhi/2); 
	
	G4LogicalVolume* logic_inner_sphere = new G4LogicalVolume(inner_sphere,
							scint_mat, "inner_sphere"); 

	G4VPhysicalVolume* physic_inner_sphere = new G4PVPlacement(nullptr,  // no rotation
						G4ThreeVector(),                     // at position
						logic_inner_sphere,              // its logical volume
						"inner_sphere",                 // its name
						logic_inner_cylinder,                 // its mother  volume
						true,                    // yes boolean operation
						0,                        // copy number
						checkOverlaps);           // overlaps checking
						
	// PMT array
	
	char pmt_name[100]="pmt";
	char pmt_logic_name[100]="pmt_logic";
	G4ThreeVector* eZ = new G4ThreeVector(0,0,1);
	G4ThreeVector* axis = new G4ThreeVector();
	int n_pmt = 20;

	G4Tubs* pmt_arr[n_pmt];
	
	double f = (1+sqrt(5))/2;
	if (n_pmt == 20){
		for (int i = 0; i < n_pmt; i++){
		double r = 800;
		double pmt_20[20][3] = {
			            {r/sqrt(3), r/sqrt(3), r/sqrt(3)},
			            {-r/sqrt(3), r/sqrt(3), r/sqrt(3)},
			            {-r/sqrt(3), -r/sqrt(3), r/sqrt(3)},
			            {-r/sqrt(3), -r/sqrt(3), -r/sqrt(3)},
			            {r/sqrt(3), -r/sqrt(3), -r/sqrt(3)},
			            {r/sqrt(3), r/sqrt(3), -r/sqrt(3)},
			            {-r/sqrt(3), r/sqrt(3), -r/sqrt(3)},
			            {r/sqrt(3), -r/sqrt(3), r/sqrt(3)},

			            {0, r/(sqrt(3)*f), r*f/sqrt(3)},
			            {0, -r/(sqrt(3)*f), r*f/sqrt(3)},
			            {0, -r/(sqrt(3)*f), -r*f/sqrt(3)},
			            {0, r/(sqrt(3)*f), -r*f/sqrt(3)},

			            {r/(sqrt(3)*f), r*f/sqrt(3), 0},
			            {-r/(sqrt(3)*f), r*f/sqrt(3), 0},
			            {-r/(sqrt(3)*f), -r*f/sqrt(3), 0},
			            {r/(sqrt(3)*f), -r*f/sqrt(3), 0},

			            {r*f/sqrt(3), 0, r/(sqrt(3)*f)},
			            {-r*f/sqrt(3), 0, r/(sqrt(3)*f)},
			            {-r*f/sqrt(3), 0, -r/(sqrt(3)*f)},
			            {r*f/sqrt(3), 0, -r/(sqrt(3)*f)}
			        };
		    //sprintf(pmt_name, "pmt_%d", i);
		    G4String* g4_pmt_name = new G4String(pmt_name);
		    
		    G4ThreeVector* pmt_vect = new G4ThreeVector(pmt_20[i][0], pmt_20[i][1], pmt_20[i][2]);
		    
		    G4ThreeVector* inversed_pmt_vect = new G4ThreeVector(-pmt_20[i][0], -pmt_20[i][1], -pmt_20[i][2]);
		    
		    G4Tubs* pmt_sample = new G4Tubs(pmt_name, 0*mm, PMTWindowRadius, PMTWindowHalfLength, StartPhi, DeltaPhi);
		    //G4cout << "angle: " << asin(101/131) << G4endl;
		 
		    
		    G4Sphere* pmt_glass = new G4Sphere("pmt_glass_sample", 110*mm, 131*mm, StartPhi, DeltaPhi, 0*deg, 50.44*deg);
		   
		    
		    G4LogicalVolume* logic_pmt_glass = new G4LogicalVolume(pmt_glass, BSG, "pmt_glass_logic");  
		    
		    //sprintf(pmt_logic_name, "pmt_%d_logic", i);
		    
		    G4LogicalVolume* logic_pmt_sample = new G4LogicalVolume(pmt_sample, outer_sphere_mat, pmt_logic_name);  
		      
		    
		    axis = vector_product(eZ, pmt_vect);
		    
		    G4double angle = angle_between(eZ, inversed_pmt_vect);
		    //G4cout << i<< " angle: "<< angle << G4endl;
		    //G4cout << i<< pmt_20[i][0] << "	" << pmt_20[i][1] << "	" << pmt_20[i][2] << G4endl;
		    
		    //if (angle > 180){ angle -= 180; }
		    G4RotationMatrix* rotmat = new G4RotationMatrix();
		    rotmat -> rotate(angle, axis); 

		    G4VPhysicalVolume* physic_pmt_sample = new G4PVPlacement(rotmat,  // no rotation
							    G4ThreeVector(pmt_20[i][0], pmt_20[i][1], pmt_20[i][2]),          // at pmt coord
							    logic_pmt_sample,                 // its logical volume
							    pmt_name,               // its name
							    logic_inner_cylinder,               // its mother  volume
							    false,                    // no boolean operation
							    0,                        // copy number
							    checkOverlaps);           // overlaps checking	
							    
		    G4VPhysicalVolume* physic_pmt_glass = new G4PVPlacement(rotmat,  // no rotation
							    G4ThreeVector(pmt_20[i][0], pmt_20[i][1], pmt_20[i][2]),          // at pmt coord
							    logic_pmt_glass,                 // its logical volume
							    "pmt_glass",               // its name
							    logic_inner_cylinder,               // its mother  volume
							    false,                    // no boolean operation
							    0,                        // copy number
							    checkOverlaps);           // overlaps checking					    	   
		//------------------------------------------------------
		// Surfaces and boundary processes
		//------------------------------------------------------

		// Scintillator - PMT window surface 
		
		G4OpticalSurface* OpScintPMTWinSurface = 
					new G4OpticalSurface("ScintPMTWinSurface");
					
		// г4 сфере умеет задваться куском чтобы сделать стекло
		
		OpScintPMTWinSurface->SetType(dielectric_dielectric);
		OpScintPMTWinSurface->SetModel(glisur);
		OpScintPMTWinSurface->SetFinish(polished);
		
		G4LogicalBorderSurface* ScintPMTWinSurface = 
	    				new G4LogicalBorderSurface("ScintPMTWinSurface",physic_inner_sphere,physic_pmt_sample,
								   OpScintPMTWinSurface);

		    pmt_arr[i] = pmt_sample;
		}
	} else { // 12 pmt
		double r = 420; // sqrt(3) factor
		double pmt_12[12][3] = {
				    {f*r,r,0},
				    {f*r,-r,0},
				    {-f*r,-r,0},
				    {-f*r,r,0},
				    
				    {r,0,f*r},
				    {-r,0,f*r},
				    {-r,0,-f*r},
				    {r,0,-f*r},
				    
				    {0,f*r,r},
				    {0,f*r,-r},
				    {0,-f*r,-r},
				    {0,-f*r,r}
				};
				
		for (int i = 0; i < n_pmt; i++){
		    //sprintf(pmt_name, "pmt_%d", i);
		    G4String* g4_pmt_name = new G4String(pmt_name);

		    G4ThreeVector* pmt_vect = new G4ThreeVector(pmt_12[i][0], pmt_12[i][1], pmt_12[i][2]);
		    
		    G4ThreeVector* inversed_pmt_vect = new G4ThreeVector(-pmt_12[i][0], -pmt_12[i][1], -pmt_12[i][2]);
		    
		    G4Tubs* pmt_sample = new G4Tubs(pmt_name, 0*mm, PMTWindowRadius, PMTWindowHalfLength, StartPhi, DeltaPhi);
		    //G4cout << "angle: " << asin(101/131) << G4endl;
		 
		    
		    G4Sphere* pmt_glass = new G4Sphere("pmt_glass_sample", 110*mm, 131*mm, StartPhi, DeltaPhi, 0*deg, 50.44*deg);
		   
		    
		    G4LogicalVolume* logic_pmt_glass = new G4LogicalVolume(pmt_glass, BSG, "pmt_glass_logic");  
		    
		    //sprintf(pmt_logic_name, "pmt_%d_logic", i);
		    
		    G4LogicalVolume* logic_pmt_sample = new G4LogicalVolume(pmt_sample, outer_sphere_mat, pmt_logic_name);  
		      
		    
		    axis = vector_product(eZ, pmt_vect);
		    
		    G4double angle = angle_between(eZ, inversed_pmt_vect);
		    //G4cout << i<< " angle: "<< angle << G4endl;
		    //G4cout << i<< pmt_12[i][0] << "	" << pmt_12[i][1] << "	" << pmt_12[i][2] << G4endl;
		    
		    //if (angle > 180){ angle -= 180; }
		    G4RotationMatrix* rotmat = new G4RotationMatrix();
		    rotmat -> rotate(angle, axis); 

		    G4VPhysicalVolume* physic_pmt_sample = new G4PVPlacement(rotmat,  // no rotation
							    G4ThreeVector(pmt_12[i][0], pmt_12[i][1], pmt_12[i][2]),          // at pmt coord
							    logic_pmt_sample,                 // its logical volume
							    pmt_name,               // its name
							    logic_inner_cylinder,               // its mother  volume
							    false,                    // no boolean operation
							    0,                        // copy number
							    checkOverlaps);           // overlaps checking	
							    
		    G4VPhysicalVolume* physic_pmt_glass = new G4PVPlacement(rotmat,  // no rotation
							    G4ThreeVector(pmt_12[i][0], pmt_12[i][1], pmt_12[i][2]),          // at pmt coord
							    logic_pmt_glass,                 // its logical volume
							    "pmt_glass",               // its name
							    logic_inner_cylinder,               // its mother  volume
							    false,                    // no boolean operation
							    0,                        // copy number
							    checkOverlaps);           // overlaps checking					    	   
		//------------------------------------------------------
		// Surfaces and boundary processes
		//------------------------------------------------------

		// Scintillator - PMT window surface 
		
		G4OpticalSurface* OpScintPMTWinSurface = 
					new G4OpticalSurface("ScintPMTWinSurface");
					
		// г4 сфере умеет задваться куском чтобы сделать стекло
		
		OpScintPMTWinSurface->SetType(dielectric_dielectric);
		OpScintPMTWinSurface->SetModel(glisur);
		OpScintPMTWinSurface->SetFinish(polished);
		
		G4LogicalBorderSurface* ScintPMTWinSurface = 
	    				new G4LogicalBorderSurface("ScintPMTWinSurface",physic_inner_sphere,physic_pmt_sample,
								   OpScintPMTWinSurface);

		    pmt_arr[i] = pmt_sample;
		}
	}
	
	//------------------------------------------------------
	// visualization attributes
	//------------------------------------------------------
	
	//logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
	
	//Green color for scintillator sphere
	G4VisAttributes* Att1= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	logic_inner_sphere->SetVisAttributes(Att1);
	
	G4VisAttributes* top_bottom_att = new G4VisAttributes(G4Colour(0.,1.,0.));
  	top_bottom_att->SetForceWireframe(true);
  	logic_bottom_lid->SetVisAttributes(top_bottom_att);
  	logic_top_lid->SetVisAttributes(top_bottom_att);
	
	G4VisAttributes* out_sphere_att = new G4VisAttributes(G4Colour(1.,1.,0.));
  	out_sphere_att->SetForceWireframe(true);
  	logic_outer_sphere->SetVisAttributes(out_sphere_att);
  	
  	G4VisAttributes* out_cylinder_att = new G4VisAttributes(G4Colour(0.,1.,1.));
  	out_cylinder_att->SetForceWireframe(true);
  	logic_outer_cylinder->SetVisAttributes(out_cylinder_att);
  	
  	logic_inner_sphere -> SetVisAttributes (false);
  	
  	logic_inner_cylinder -> SetVisAttributes (false);
	
	
    
	//
	// always return the physical World
	//
	
	return physiWorld;
}

}
