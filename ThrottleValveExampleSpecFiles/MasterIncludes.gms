*****************************
***** MasterInclude.gms *****
*****************************

* What to specify different include files? Edit this file!

* Generic flowsheet parameters
$setglobal GenericParamsFile "ThrottleValveExampleSpecFiles/GenericParameters.gms";

* Flowsheet topology
$setglobal FlowsheetTopologyFile "ThrottleValveExampleSpecFiles/FlowsheetTopology_example2.gms";

* Thermodynamic Data
$setglobal ThermoDataFile "ThrottleValveExampleSpecFiles/Thermo_Data.gms";

* Problem specific bounds
$setglobal ProbSpecificBounds "ThrottleValveExampleSpecFiles/BoundsExampleOnly.gms";

* Problem specific objective function(s)
$setglobal ProbObjFncsFile "ThrottleValveExampleSpecFiles/ProbSpecObjFunctions.gms";

* Initialize problem specific objective function(s)
$setglobal InitProbObjFncsFile "ThrottleValveExampleSpecFiles/ProbObjSpecInit.gms";

* Default values for multistart initialization
$setglobal ExternalInitDefaultFile "ThrottleValveExampleSpecFiles/ExternalInitValues.gms";

* Declarations for Problem/Section Specific Initialization
$setglobal ProbSecSpecDeclFile "ThrottleValveExampleSpecFiles/ProbSecSpecDecl.gms";

* Routines for Problem/Section Specific Initialization
$setglobal ProbSecSpecInitFile "ThrottleValveExampleSpecFiles/ProbSecSpecInit.gms";

* Routines for Problem/Objective Function Specific Initialization
$setglobal ProbObjSpecInitFile "ThrottleValveExampleSpecFiles/ProbObjSpecInit.gms";

* RDX file containing flowrates, pressures and temperatures.
* Used for initialization
$setglobal SolForInit "'FlowsheetOpt_TrayByTrayb_SavedSol18'";

* Problem specific pruning
$setglobal ProbSpecPruneFile "ThrottleValveExampleSpecFiles/ProbSpecPruning.gms";
