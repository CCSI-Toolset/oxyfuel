**********************************
***** ExternalInitValues.gms *****
**********************************

* This file contains values for all of the external multistart initialization
* parameters. These values overwrite the imported values when "ExternalInit"
* is set to 0.

if(ExternalInit eq 0,
  NumStageInit = 25;
  PhiLBCEOS = -8;
  UBASFact = 10;
  LBASFact = 0.0001;
  Tmin = 65;
  Prune = 0;
  alpha = -0.1;
  HRATSimple = 6;
  ProbSpecParam1 = 0;
  ProbSpecParam2 = 0.55;

  HRATPlusBD = 0;

  FirstSolverSimple = 1;
  FirstSolverCEOS = 1;
  FirstSolverTray = 1;

  HeatIntegModel = 1;

  FirstDerivEpsilon = -7;
  Inter1Epsilon = -7;
  CascInterLB = -7;
  PlusStages = 10;
  FixEffForInit = 0;
  SelectCEOSModel = 3;
  PruneConfig = 0;

  ProbSpecParam3 = 0.95;

  HRATCEOS = 1.5;

  SimpleCscInitProb = 1;
  LoadInitPoint = 1;
);

Parameters
  O2RecLowDelta,
  O2RecLowCEOS,
  o2pureSpec;

O2RecLowDelta = ProbSpecParam1;
O2RecLowCEOS = ProbSpecParam2;
o2pureSpec = ProbSpecParam3;

