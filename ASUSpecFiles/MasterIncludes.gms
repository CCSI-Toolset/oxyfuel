*****************************
***** MasterInclude.gms *****
*****************************

* What to specify different include files? Edit this file!

* Generic flowsheet parameters
$setglobal GenericParamsFile "ASUSpecFiles/GenericParameters.gms";

* Flowsheet topology
$setglobal FlowsheetTopologyFile "ASUSpecFiles/FlowsheetTopology_ASU_NoDualReboiler.gms";

* Thermodynamic Data
$setglobal ThermoDataFile "ASUSpecFiles/Thermo_Data.gms";

* Problem specific bounds
$setglobal ProbSpecificBounds "ASUSpecFiles/BoundsASUOnly.gms";

* Problem specific objective function(s)
$setglobal ProbObjFncsFile "ASUSpecFiles/ProbSpecObjFunctions.gms";

* Initialize problem specific objective function(s)
$setglobal InitProbObjFncsFile "ASUSpecFiles/ProbObjSpecInit.gms";

* Default values for multistart initialization
$setglobal ExternalInitDefaultFile "ASUSpecFiles/ExternalInitValues.gms";

* Declarations for Problem/Section Specific Initialization
$setglobal ProbSecSpecDeclFile "ASUSpecFiles/ProbSecSpecDecl.gms";

* Routines for Problem/Section Specific Initialization
$setglobal ProbSecSpecInitFile "ASUSpecFiles/ProbSecSpecInit.gms";

* Routines for Problem/Objective Function Specific Initialization
$setglobal ProbObjSpecInitFile "ASUSpecFiles/ProbObjSpecInit.gms";

* RDX file containing flowrates, pressures and temperatures.
* Used for initialization
$setglobal SolForInit "'FlowsheetOpt_TrayByTrayb_SavedSol18'";

* Problem specific pruning
$setglobal ProbSpecPruneFile "ASUSpecFiles/ProbSpecPruning.gms";
