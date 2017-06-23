* TOP OF FILE

* This is the "main" file for simple LPC synthesis and optimization. Include files are
* extensively used.
*
* Note: This file is based on ASU27_DualReboiler
*
* The optimization problem is solved in six stages/sections.
*
* 0. Either load stream data from a file or initial with constant stream properties
*
* 1. Simple thermodynamic flowsheet optimziation
*    In the first stage a simplified flowsheet model featuring correlation based
*    thermodynamic equations is optimized and heat integrated. This determines
*    reasonable temperatures, pressures, flows, etc for all of the ASU streams.                                                      )
*
* 2. Cubic CEOS model initialization
*    In the second stage stream properties are fixed, and the CEOS model
*    are initialized. This is done by solving three optimization problems.
*    During this stage the shadow streams used for bubble/dew point calculations
*    are "activated".
*
* 3. CEOS flowsheet model optimization with Edmister cascades
*    In the third stage the ASU with CEOS thermo is optimized and heat
*    integrated. The group method Edmister-based cascade approximation is used
*    in this model, which results in the number of theoretical trays being a
*    continous optimization variable. The optimization parameter is originally
*    solved with a large smoothing epsilon and loose group method bounds.
*    Epsilon is decrease and the bounds are tightened in subsequent solves. In
*    the final solves the number of stages in the cascades is rounded and fixed.
*
* 4. CEOS flowsheet model optimization with tray-by-tray cascades
*    To facilitate model validation, the Edmister cascade model is replaced with
*    a tray-by-tray (MESH) model in the fourth stage. For these final solved the number
*    of trays remains fixed. Conversion from the group method to cascade model
*    is automated with dynamic sets in GAMS. Initialization is also automated
*    with loops.
*
* 5. The constant heat capacity assumption is checked by decomposing heat exchangers
*    into multiple subunits
*
* 6. CEOS flowsheet model is re-optimized with MESH cascade model and refined heat exchange units

$setglobal GlobalIterlim "10000";

* global options
$eolcom //
option optcr = 1E-10; ;
*option solprint = ON ;
*option reslim = 600;
option reslim = 3600;
option iterlim = %GlobalIterlim%;
$onmacro
$onempty
$onexpand

$set matout "'gamsstat.gdx', Sec1Stat, Sec2Stat, Sec3Stat, Sec4Stat, Sec5Stat, Sec6Stat ";
// end of global options

***** Stage Mode Settings *****
* These parameters control the "operational mode" for each stage.

Parameter
ProblemSimpleFlowsheetOpt         Initialization mode                                     /1/
* Mode 0 - skip. don't solve
* Mode 1 - solve

ProblemInitCEOS          Initialization mode                                     /1/
* Mode 0 - don't solve
* Mode 1 - solve
* Mode 2 - don't solve and load previous solution

ProblemCEOSFlowsheetOpt           Initialization mode                                     /1/
* Mode 0 - don't solve
* Mode 1 - solve
* Mode 2 - don't solve and load previous solution

ProblemTrayByTray        Solve mode                                              /0/
* Mode 0 - don't solve
* Mode 1 - solve
* Mode 2 - don't solve and load previous solution

ProblemHeatChecker              Solve mode                                              /0/
* Mode 0 - don't solve
* Mode 1 - solve
* Mode 2 - don't solve and load previous solution

ProblemTrayByTrayDecomp  Solve mode                                              /0/
* Mode 0 - don't solve
* Mode 1 - solve

ThermoSwitch       Toggle between thermo methods                           /1/
* Mode 1 - simple thermo (vapor pressure correlations)
* Mode 2 - CEOS thermo
* Should start as 1... don't change here

SectionSwitch      Toggle between sections for problem specific init.            /0/;
* Should start at 0... don't change here


***** Specify Master Include File *****
$include ../ASUSpecFiles/MasterIncludes.gms
*$include ../NewSpecFiles/MasterIncludes.gms
$include ../ThrottleValveExampleSpecFiles/MasterIncludes.gms

***********************************************************************
***** Do not change anything below this line for normal operation *****
***********************************************************************

**************************************
***** Setup Solve Status Storage *****
**************************************

Set stat /modelstat,solvestat,timer/;
Parameters
  Sec1Stat(stat)
  Sec2Stat(stat)
  Sec3Stat(stat)
  Sec4Stat(stat)
  Sec5Stat(stat)
  Sec6Stat(stat);

***********************************
***** External Initialization *****
***********************************

Parameters
  NumStageInit   Number of initial stages in cascades
  PhiLBCEOS      Lower bound for phi in cascade model for CEOS problems
  UBASFact       Upper bound on absorption or stripping factor
  LBASFact       Lower bound on absorption or stripping factor
  Tmin           Min temperature (K) for the flowsheet
  Prune          Manually remove equipment before FlowsheetOpt_Simple?
  alpha          Small number used for heat integration phase changes
  saveAlpha
  HRATSimple     HRAT for FlowsheetOpt_Simple
  HRATCEOS       HRAT for FlowsheetOpt_CEOS and TrayByTray
  HRATPlusBD     Plus value for HRAT before decomposition of heat exchange units
  FirstSolverSimple
  FirstSolverCEOS
  FirstSolverTray
  HeatIntegModel
  FirstDerivEpsilon
  Inter1Epsilon
  CascInterLB
  PlusStages
  FixEffForInit
  SelectCEOSModel
  ExternalInit  ExternalInitialization flag                           /0/
  LinuxOS       Running on Unix or Linux?                             /0/
  LoadInitPoint Load initial flows-temps-pressures from file
  ProbSpecParam1
  ProbSpecParam2
  ProbSpecParam3
  SimpleCscInitProb
  CEOSCscInitProb                                                     /0/
  PruneConfig
  DualReboilerASU                                                     /0/
;

$if %system.filesys% == UNIX $set LinuxOS 1
display LinuxOS;

$gdxin initData
$load NumStageInit
$load PhiLBCEOS
$load UBASFact
$load LBASFact
$load Tmin
$load Prune
$load alpha
$load HRATSimple
*$load O2RecLowDelta
*$load O2RecLowCEOS
$load HRATPlusBD
$load HRATCEOS
$load FirstSolverSimple
$load FirstSolverCEOS
$load FirstSolverTray
$load HeatIntegModel
$load FirstDerivEpsilon
$load Inter1Epsilon
$load CascInterLB
$load PlusStages
$load FixEffForInit
$load SelectCEOSModel
$load PruneConfig
*$load o2pureSpec
$load SimpleCscInitProb
$load LoadInitPoint

$load ProbSpecParam1
$load ProbSpecParam2
$load ProbSpecParam3
$gdxin

$include ../%ExternalInitDefaultFile%


************************************
***** Include flowsheet models *****
************************************

* Generic parameters, some of which are used in the models
$include ../%GenericParamsFile%

* Flowsheet Topology
$include ../%FlowsheetTopologyFile%
$include ../ProcsFiles/FlowsheetTopologyProcessing.gms

* Define smoothed max operator
$macro smmax(x) ( 0.5*sqrt(Power(x,2) + epsiSmax) + 0.5*(x) )

* Load thermodynamic data and constants
$include ../%ThermoDataFile%

* Stream model
$include ../ModelFiles/StreamModel.gms

* Equipment models independent of thermo
$include ../ModelFiles/ThrmEModel.gms
$include ../ModelFiles/FlashModel.gms
$include ../ModelFiles/HeatExchangerModel.gms
$include ../ModelFiles/ReboilerModel.gms
$include ../ModelFiles/ValveModel.gms
$include ../ModelFiles/MixerModel.gms

* Thermo models
$include ../ModelFiles/CommonThermoModel.gms
$include ../ModelFiles/SimpleThermoModel.gms
$include ../ModelFiles/CEOSThermoModel.gms
$include ../ModelFiles/CEOSDewBubModel2.gms

* Equipment models dependent on thermo
$include ../ModelFiles/CondenserModel.gms
$include ../ModelFiles/CascadeModel.gms
$include ../ModelFiles/EntropyBalModel.gms
$include ../ModelFiles/SplitterModel.gms

* Heat integration model
$include ../ModelFiles/HeatIntegrationModelCommon.gms
$include ../ModelFiles/HeatIntegrationModelUpdated.gms
*$include ../ModelFiles/HeatIntegrationModelRevised.gms

* Advanced distillation models
$include ../ModelFiles/TrayByTrayModel.gms
$include ../ModelFiles/MESHwBypassModel3.gms

*--------------------------------*
*     Dummy variables
*--------------------------------*
Variables
Z ;

$include ../%ProbSecSpecDeclFile%
$include ../ModelFiles/ObjFncsGeneric.gms
$include ../%ProbObjFncsFile%

*0000000000000000000000000000000000*
************************************
********** BEGIN SECTION 0 *********
************************************
*0000000000000000000000000000000000*
SectionSwitch = 0;

Parameters
  loopExp /1/
  loopCount /1/;

*--------------------------------*
*  Load initial point from file
*--------------------------------*

F.l(Str) = 0;
Tscaled.l(Str) = 0;
P.l(Str) = 0;

*$ontext
$ontext
$gdxin %SolForInit%
$loadr F
$loadr Fc
$loadr Xc
$loadr P
$loadr Tscaled
$gdxin
$offtext

execute_load %SolForInit% F,Fc,Xc,P,Tscaled;

F.m(Str) = 0;
Fc.m(Str,Comp) = 0;
Xc.m(Str,Comp) = 0;
Tscaled.m(Str) = 0;

F.lo(Str) = -Inf;
F.up(Str) = Inf;

Fc.lo(Str,Comp) = -Inf;
Fc.up(Str,Comp) = Inf;

P.lo(Str) = -Inf;
P.up(Str) = Inf;

Tscaled.lo(Str) = -Inf;
Tscaled.up(Str) = Inf;

Tscaled.l(Str)$(Tscaled.l(Str) eq 0) = FeedT/Tref;
P.l(Str)$(P.l(Str) eq 0) = FeedP;
*F.l(Str) = max(F.l(Str), 1E-6);

T.l(Str) = Tscaled.l(Str)*Tref;

Xc.l(Str,Comp)$(SUM(Comp2, Xc.l(Str,Comp2)) < 0.1) = 1/card(Comp);

*$offtext

$ontext
F.l(Str) = 0;
Tscaled.l(Str) = 0;
P.l(Str) = 0;
$offtext

Parameters FeedTTmp /100/;

* Initialize new streams (not in the initPoint file)
loop(Str$(F.l(Str) eq 0 AND Tscaled.l(Str) eq 0 AND P.l(Str) eq 0),
  Fc.l(Str,Comp) = FeedFlow(Comp);
  F.l(Str) = SUM(Comp, Fc.l(Str,Comp));
  Xc.l(Str,Comp) = Fc.l(Str,Comp)/Max(F.l(Str),epsi);
  Tscaled.l(Str) = FeedTTmp/Tref;
  P.l(Str) = FeedP;
);


$include ../%ProbSecSpecInitFile%

EpsilonZ = FirstDerivEpsilon;
EpsilonA = Inter1Epsilon;
Inter1.lo(Str) = 10**EpsilonA;

$include ../ModelFiles/SetupPostProcessing.gms

$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%


*------------------------**
* Initialize Flows
*------------------------**

***** Option 1: Load From File *****
* In this option temperatures, pressures, compositions and flowrates
* are loaded from a file. This is analogous to starting from a known solution,
* such as results obtained from a commerical flowsheet simulator.

* Do nothing. This was already done above.

***** Option 2: Load a general guess ****
* Using this option the flowsheet is initialized to contain non-zero flows with
* some guess temperatures and pressures. In the current state, option 1 tends to
* be more effective. Refinement of the guesses in option 2 should improve
* performance.


Parameter FeedXc(Comp);

FeedXc(Comp) = FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));
Display FeedXc;

Equation EqInitObj;

EqInitObj..
  Z =e= Sum(Str, (F(Str) - 1)*(F(Str) - 1)
         + 10*SUM(Comp, (Xc(Str,Comp) - FeedXc(Comp))*(Xc(Str,Comp) - FeedXc(Comp))));


Model MassBalEqns /EqTotMolStr, EqSplitFrac, EqTotCompMolBalEd1, EqCompMolBalGnrlE,
                  EqMolBalSptr, EqMolFracSptr, PurityRecoveryEquations/;

Model MassBal /MassBalEqns, EqInitObj,EqSumFrac/;

if(ProblemSimpleFlowsheetOpt eq 1,

*  MassBal.savepoint = 1;
  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = GAMSCHK;
  );

* EITHER: 1. Solve mass balance problem to initialize most streams with non-zero flows
* OR: 2. Don't overwrite initial point (F, Fc, Tscaled, P) loaded from file
  if(LoadInitPoint eq 0,

    Fc.l(Str,Comp) = FeedFlow(Comp);
    F.l(Str) = SUM(Comp, Fc.l(Str,Comp));
    Xc.l(Str,Comp) = Fc.l(Str,Comp)/Max(F.l(Str),epsi);

    F.up(ExtraClmnFeeds) = 0;

    Solve MassBal using NLP minimizing Z;

    T.l(Str) = FeedTTmp;
    P.l(Str) = FeedP;
    T.l(FeedStr) = 300;
  );
);

$include ../%ProbSpecPruneFile%

F.fx(ManualRemove) = 0;
Fc.fx(ManualRemove,Comp) = 0;

F.l(Str)$(NOT ManualRemove(Str)) = max(F.l(Str), 1E-3);

$include ../ProcsFiles/ManualRemoveEquipment.gms

if(DebugMode > 0,
  display ManualRemove;
);

* Idea for initializing equipment that are "cut out" from the loaded point
* Written May 15th, 2014

$ontext
loop(GnrlE$(NOT InactiveGnrlE(GnrlE)),
  if( SUM(Str$(NOT InactiveStr(Str) AND InGnrlE(GnrlE,Str)), F.l(Str)) < epsi*epsi,
    loop(Str$(NOT InactiveStr(Str) AND InGnrlE(GnrlE,Str)),
      F.l(Str) = 0.01;
      Xc.l(Str,Comp) = FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));
      Fc.l(Str,Comp) = F.l(Str)*Xc.l(Str,Comp);
    );
  );
);
$offtext

scalar timer;

*1111111111111111111111111111111111*
************************************
********** BEGIN SECTION 1 *********
************************************
*1111111111111111111111111111111111*

SectionSwitch = 1;
timer = timeElapsed;

$include ../%ProbSecSpecInitFile%

*********************************************
***** General Initialization and Bounds *****
*********************************************

***** General Specs *****

Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;
Tscaled.l(Str) = T.l(Str)/Tref;

* Only execute these scripts if P, Tscaled, F, etc not loaded from a file
if(LoadInitPoint eq 0,
* Initialize zones - depreciated
* $include ../IncludeFiles/InitZones.gms

* Set heat exchanger temperatures and pressures - requires debugging
* $include ../IncludeFiles/InitHtEx.gms
);

***** Initialize Thermo *****

$include ../InitFiles/InitSmpThermo.gms

***** Initialize Flash *****

Qin.fx(Flsh) = 0;
Qout.fx(Flsh) = 0;

PFlsh.l(Flsh) = SUM(Str$OutLGnrlE(Flsh,Str), P.l(Str));
TFlsh.l(Flsh) = SUM(Str$OutLGnrlE(Flsh,Str), Tscaled.l(Str));

***** Initialize Valves *****

Qin.fx(Valve) = 0;
Qout.fx(Valve) = 0;

$include ../InitFiles/Init_Smp_Valve.gms

***** Initialize Partial Reboilers *****

Qout.fx(PReb) = 0;

PPReb.l(PReb) = SUM(Str$OutLGnrlE(PReb,Str), P.l(Str));

***** Initialize Total Condenser *****

PTCond.l(TCond) = SUM(Str$OutTCond(TCond,Str), P.l(Str));
PTCond.up(TCond) = Pmax;

***** Initialize Splitters *****
$include ../InitFiles/InitSplitter.gms

***** Initialize Q in GnrlE *****

$include ../InitFiles/InitGnrlE.gms

***** Initialize Cascade *****

VNEd1.l(Ed1) = SUM(Str $InVEd1(Ed1,Str), F.l(Str));
L1Ed1.l(Ed1) =  SUM(Str $InLEd1(Ed1,Str), F.l(Str)) ;

* Stages per cascade
StgEd1.l(Ed1) = NumStageInit ;

$include ../InitFiles/Init_Cascade.gms

***** Initialize Objective *****

Model FlowsheetOpt /EqStreams,EqCascade,EqThrmE,EqFlsh, EqPReb,EqValve,EqTCond,EqSpltr,EqObj/;

Model FlowsheetOpt_Simple /FlowsheetOpt, SimpleThermo, EqSpltrSimpleThermo/;

*Model FlowsheetOpt_Simple6_R /FlowsheetOpt_Simple, EqCsdSimple, EqHtEx, HeatIntegRevised, obj9/;
Model FlowsheetOpt_Simple6_U /FlowsheetOpt_Simple - EqTscaled, EqCsdSimple, EqHtEx, HeatIntegUpdated, ObjSec1/;

**************************************************
***** Special Cascade Initialization Problem *****
**************************************************

$include ../ModelFiles/InitializationEquations.gms

count = 0;

if(ProblemSimpleFlowsheetOpt eq 1 AND SimpleCscInitProb eq 1 AND SUM(Ed1$(NOT InactiveCsc(Ed1)), 1) > 0,

  Fc.l(CscStr,Comp) = F.l(CscStr)*Xc.l(CscStr,Comp);

* Setup InitEqns for InitMassBal and InitCscNew

  InitStr1(Str) = no;
  InitStr1(CscStr) = yes;

  InitStr2(Str) = no;
  InitStr2(Str)$(NOT CscStr(Str)) = yes;
$include ../InitFiles/InitInitEqns.gms

  option NLP = GAMSCHK;

  solve InitMassBal using NLP minimizing FcPen;

$include ../InitFiles/InitSmpThermo.gms

  InitCscNewObj.iterlim = 1500;
  solve InitCscNewObj using NLP minimizing Z;

  if(ProblemInitCEOS eq 0,
    option NLP = CONVERTD;
    InitCscNewObj.optfile = 1;
    solve InitCscNewObj using NLP minimizing Z;
  );
);

ComplPen = 10;

***** Initialize Heat Integration *****

if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
else
*$include ../InitFiles/Init_HtInteg_Revised.gms
);

sL.fx(LiqStr) = 0;
sV.fx(VapStr) = 0;

if(ProblemSimpleFlowsheetOpt eq 1,

***** Initialize Slacks *****
  liqPen.l = sum(ThrmE, sum(Str$OutLGnrlE(ThrmE,Str), sL.l(Str)*F.l(Str)));
  vapPen.l = sum(ThrmE, sum(Str$OutVGnrlE(ThrmE,Str), sV.l(Str)*F.l(Str)));

  ComplPen = 10;
  epsiSmax = 1E-4;


$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%

$include ../InitFiles/ThrmEBubDew.gms
$include ../%ProbObjSpecInitFile%

  ComplPen = 1;

  if(FirstSolverSimple eq 1,
    option NLP = GAMSCHK;
  elseif FirstSolverSimple eq 2,
    option NLP = CONOPT;
*    FlowsheetOpt_Simple6_R.optfile = 1;
    FlowsheetOpt_Simple6_U.optfile = 1;
  elseif FirstSolverSimple eq 3,
    option NLP = IPOPTH;
  elseif FirstSolverSimple eq 4,
    option NLP = MINOS;
  elseif FirstSolverSimple eq 5,
    option NLP = SNOPT;
  else
    option NLP = CONOPT;
  );

*  FlowsheetOpt_Simple6_R.savepoint=1;
*  FlowsheetOpt_Simple6_U.savepoint=1;

*  FlowsheetOpt_Simple6_U.iterlim = 30000;

epsiSmax = 1E-3;

*  Tscaled.l(Str) = T.l(Str)/Tref;
*  Tscaled.lo(Str) = T.lo(Str)/Tref;
*  Tscaled.up(Str) = T.up(Str)/Tref;

*  FlowsheetOpt_Simple6_U.reslim = 3*60;

for(loopCount = 1 to 1 by 1,
  for(loopExp = -1 to 0 by 1,

    ComplPen = Power(10, loopExp+1+StartCPPower);

    CascInterA.lo(Ed1,Comp) = Power(10, CascInterLB - 1 - loopExp);
    CascInterS.lo(Ed1,Comp) = Power(10, CascInterLB - 1 - loopExp);

    epsiSmax = Power(10, -5 - loopExp);

* Initialize Heat Integration
    if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
    else
*$include ../InitFiles/Init_HtInteg_Revised.gms
    );

    if(HeatIntegModel eq 1,
      solve FlowsheetOpt_Simple6_U using NLP minimizing Z9;
    else
*     solve FlowsheetOpt_Simple6_R using NLP minimizing Z9;
    );

    T.l(Str) = Tscaled.l(Str)*Tref;
    display T.l, ComplPen, epsiSmax;

    ChckCmpl(0);
*$include ../ProcsFiles/RemoveUnusedEquipment.gms

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = GAMSCHK;
    );
*    FlowsheetOpt_Simple6_R.optfile = 0;
    FlowsheetOpt_Simple6_U.optfile = 0;
  );
);

  if(HeatIntegModel eq 1,
    Sec1Stat('modelstat') = FlowsheetOpt_Simple6_U.modelstat;
    Sec1Stat('solvestat') = FlowsheetOpt_Simple6_U.solvestat;
  else
*  SimpleStat('modelstat') = FlowsheetOpt_Simple6_R.modelstat;
*  SimpleStat('solvestat') = FlowsheetOpt_Simple6_R.solvestat;
  );

  if(ProblemInitCEOS eq 0,
    option NLP = convertd;
*    FlowsheetOpt_Simple6_R.optfile = 1;
    FlowsheetOpt_Simple6_U.optfile = 1;
*    FlowsheetOpt_Simple6_R.savepoint = 0;
    FlowsheetOpt_Simple6_U.savepoint = 0;

    if(HeatIntegModel eq 1,
      solve FlowsheetOpt_Simple6_U using NLP minimizing Z9;
    else
*      solve FlowsheetOpt_Simple6_R using NLP minimizing Z9;
    );
  );

);

if( ProblemSimpleFlowsheetOpt eq 2,
  F.l(Str) = 0;
  Fc.l(Str,Comp) = 0;
  execute_loadpoint 'Sec1_SimpleFlowsheetResults.gdx';
  epsiSmax = 1E-6;
);

*$include ../ProcsFiles/AddUnusedEquipment.gms

  F.fx(ManualRemove) = 0;

* The following line dramatically impacts the solution
  F.l(Str)$(NOT ManualRemove(Str)) = max(F.l(Str), 10E-5);

* Tried this with initPoint4. No change when uncommented
* $include ../ProcsFiles/RemoveUnusedEquipment.gms

if(ProblemSimpleFlowsheetOpt eq 1,
  execute_unload 'Sec1_SimpleFlowsheetResults.gdx';
);

Sec1Stat('timer') = timeElapsed - timer;
timer = timeElapsed;

*2222222222222222222222222222222222*
************************************
********** BEGIN SECTION 2 *********
************************************
*2222222222222222222222222222222222*

SectionSwitch = 2;
$include ../%ProbSecSpecInitFile%

* Switch from Simple Thermo to CEOS Thermo Model
ThermoSwitch = 2;

***** Initialize shadow stream compositions *****

Parameters
  dPvapdT(Str,Comp)
  dTs(Str)
  Pbtemp(Str)
  Pdtemp(Str);

Scalar iter;

Set
  PrmShd(Str)   Shadow stream that is the same phase as actual process stream it's paired with / /;

* Define parameters for analytical ZEOS initialization
$include ../ProcsFiles/SetupInitZEOS.gms

if(ProblemInitCEOS eq 1,

* Ensure shadow streams are properly setup
$include ../ProcsFiles/ProcessShadowStreams.gms

$include ../InitFiles/InitializeShadowStreams.gms

* Initialize ZEOS using analytical solutions
$batinclude ../InitFiles/Init_CEOS_Basic4.gms "(fEOSStr(Str) AND ((RealStr(Str) AND NOT InactiveStr(Str)) OR ActShdStr(Str)))"

* Set Xc, T and P targets for objective function
InitStr1(Str) = no;
InitStr1(Str)$(fEOSStr(Str) AND ((RealStr(Str) AND NOT InactiveStr(Str)) OR ActShdStr(Str))) = yes;

InitStr2(Str) = no;

$include ../InitFiles/InitInitEqns.gms

MinSlacks(Str) = no;

MinSlacks(Str)$(RealStr(Str) AND fEOSStr(Str) AND NOT InactiveStr(Str) AND F.l(Str) > 0) = yes;
MinSlacks(Str)$(fEOSStr(Str) AND ActShdStr(Str)) = yes;

display MinSlacks, Ttarget;


Z.l = SUM(MinSlacks, sV.l(MinSlacks) + sL.l(MinSlacks)) + 100*TPen.l + 10*XcPen.l + 10*PPen.l;

***** Store initialization *****

InitCEOS5.iterlim = 0;
*InitCEOS5.savepoint = 1;

option NLP = GAMSCHK;
Solve InitCEOS5 using NLP minimizing Z;

***** Try with CONOPT *****

* Initialize basic CEOS equations (everything required for calculating ZEOS)
InitCEOS5.iterlim = %GlobalIterlim%;
InitCEOS5.savepoint = 0;
option NLP = GAMSCHK;

if(ProblemInitCEOS eq 1,
  Solve InitCEOS5 using NLP minimizing Z;
else
  InitCEOS5.modelstat = -1;
);
***** Use IPOPT if CONOPT fails *****
if(ProblemInitCEOS eq 1 AND InitCEOS5.modelstat > %ModelStat.Infeasible% - 1 AND InitCEOS5.modelstat < %ModelStat.Intermediate Infeasible% + 1,
  execute_loadpoint 'InitCEOS5_p';
  option NLP = IPOPTH;
  Solve InitCEOS5 using NLP minimizing Z;
);

InitStr1(Str) = no;
InitStr1(Str)$(fEOSStr(Str) AND (RealStr(Str) AND NOT InactiveStr(Str))) = yes;
$include ../InitFiles/InitInitEqns.gms

* Initialize departure functions
$include ../InitFiles/Init_CEOS_Adv.gms

* Solve with departure functions and phase stability, dew/bubble point constraints
  option NLP = GAMSCHK;

*option iterlim = 0;
*Solve InitCEOS6 using NLP minimizing Z;

  option iterlim = %GlobalIterlim%;
  Solve InitCEOS6 using NLP minimizing Z;

  if(Z.l > 10,
    option NLP = KNITRO;
    Solve InitCEOS6 using NLP minimizing Z;
  );

  Sec2Stat('modelstat') = InitCEOS6.modelstat;
  Sec2Stat('solvestat') = InitCEOS6.solvestat;

if(ProblemCEOSFlowsheetOpt eq 0,
  option NLP = CONVERTD;
  InitCEOS6.optfile = 1;
  Solve InitCEOS6 using NLP minimizing Z;
);


);

*if(ProblemInitCEOS eq 2,
*  execute_loadpoint 'InitCEOS6_p.gdx';
*);

if(ProblemInitCEOS eq 1,
  Z9.l = Z.l;
  execute_unload 'Sec2_InitCEOSResults.gdx';

elseif ProblemInitCEOS eq 2,
  execute_load 'Sec2_InitCEOSResults.gdx';
);

Sec2Stat('timer') = timeElapsed - timer;
timer = timeElapsed;

*3333333333333333333333333333333333*
************************************
********** BEGIN SECTION 3 *********
************************************
*3333333333333333333333333333333333*

SectionSwitch = 3;

* Adjust zone approach temperatures
  HRAT(HIZone) = HRATCEOS + HRATPlusBD;

* This line is neccissary with the JT specific valve model, but has no effect
* with the H(T,P) simple thermo correlation
CalcGnrlE(Valve) = yes;

Model FlowsheetOpt_CEOS_Base /FlowsheetOpt - EqTscaled, CubicEOSEqns, EqCsdCEOS, BubbleDewEqns, PhaseStabilityEqns, EqHtEx, EqSpltrCEOSThermo/;

*Model FlowsheetOpt_CEOS9b2_R /FlowsheetOpt_CEOS_Base, HeatIntegRevised, obj9, /;
Model FlowsheetOpt_CEOS9b2_U /FlowsheetOpt_CEOS_Base, HeatIntegUpdated, ObjSec3/;

Equation EqInitSimpleThermoObj;

*********************************
*** Initialization and Bounds ***
*********************************

***** General *****

*PmaxC('HPC') = 10;
*PmaxC('LPC') = 5;

$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%

$include ../%ProbSecSpecInitFile%


Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;

F.fx(ManualRemove) = 0;

* Init HIG
$include ../InitFiles/Init_CEOS_Adv.gms

***** Cascade *****

$include ../InitFiles/Init_Cascade.gms


***** Special Cascade Initialization Problem *****

* Setup InitEqns for InitCscCEOS

  InitStr1(Str) = no;
  InitStr1(CscStr) = yes;

  InitStr2(Str) = no;
  InitStr2(Str)$(NOT CscStr(Str)) = yes;
$include ../InitFiles/InitInitEqns.gms

epsiSmax = 1E-3;

if(CEOSCscInitProb eq 1,

  loopExp = -2;
  if(PhiLBCEOS eq 0,
    PhiAEd1.lo(Ed1,Comp) = 0;
    PhiSEd1.lo(Ed1,Comp) = 0;
  else
    PhiAEd1.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);
    PhiSEd1.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);

    DummyAe.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);
    DummySe.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);
  );
  CascInterA.lo(Ed1,Comp) = Power(10, CascInterLB - loopExp);
  CascInterS.lo(Ed1,Comp) = Power(10, CascInterLB - loopExp);
  EpsilonZ = Power(10, FirstDerivEpsilon - loopExp);
  Inter1.lo(Str) = Power(10, Inter1Epsilon - loopExp);

  if(ProblemCEOSFlowsheetOpt > 0 AND ProblemCEOSFlowsheetOpt < 3,
    option NLP = GAMSCHK;
*  option NLP = IPOPTH;
*    InitCscCEOS.reslim = 45;
    InitCscCEOS.iterlim = 1500;
    solve InitCscCEOS using NLP minimizing Z;

  );

  StgEd1.lo(Ed1) = LBStgEd1 ;
  StgEd1.up(Ed1) = UBStgEd1 ;
);

$include ../InitFiles/InitGnrlE.gms

if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
else
*$include ../IncludeFiles/Init_HtInteg_Revised.gms
);

$include ../%InitProbObjFncsFile%
$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%

if(ProblemCEOSFlowsheetOpt eq 1,

  ComplPen = 1;

$ontext
  CmprPwr.l(FeedStr) = F.l(FeedStr)*Ncmpr*gamma/(gamma - 1)*Rsi*300*
         ( P.l(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1)
         /(Fc.l('S15','O2') + Fc.l('S17','O2'))/(1000*115.2);
$offtext

  PhiAEd1.lo(Ed1,Comp) = Power(10,PhiLBCEOS);
  PhiSEd1.lo(Ed1,Comp) = Power(10,PhiLBCEOS);

$include ../InitFiles/ThrmEBubDew.gms

if(Prune eq 1,
$include ../ProcsFiles/RemoveUnusedEquipment.gms
);

*  FlowsheetOpt_CEOS9b2_R.savepoint = 1;
*  FlowsheetOpt_CEOS9b2_U.savepoint = 1;

  if(FirstSolverCEOS eq 1,
    if(ExternalInit eq 1,
      option NLP = CONOPT;
      FlowsheetOpt_CEOS9b2_U.optfile = 0;
    else
      option NLP = GAMSCHK;
    );
  elseif FirstSolverCEOS eq 2,
    option NLP = CONOPT;
    FlowsheetOpt_CEOS9b2_U.optfile = 1;
*    FlowsheetOpt_CEOS9b2_R.optfile = 1;
  else
    option NLP = IPOPTH;
  );

epsiSmax = 1E-3;

*  Tscaled.l(Str) = T.l(Str)/Tref;
*  Tscaled.lo(Str) = T.lo(Str)/Tref;
*  Tscaled.up(Str) = T.up(Str)/Tref;

*option iterlim = 0;

for(loopCount = 1 to 1 by 1,
*  for(loopExp = -1 to 1 by 1,
  for(loopExp = -1 to -1 by 1,
    ComplPen = Power(10, loopExp+1+StartCPPower);
    if(PhiLBCEOS eq 0,
      PhiAEd1.lo(Ed1,Comp) = 0;
      PhiSEd1.lo(Ed1,Comp) = 0;
    else
      PhiAEd1.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);
      PhiSEd1.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);

      DummyAe.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);
      DummySe.lo(Ed1,Comp) = Power(10, PhiLBCEOS - loopExp);

    );
    CascInterA.lo(Ed1,Comp) = Power(10, CascInterLB - loopExp);
    CascInterS.lo(Ed1,Comp) = Power(10, CascInterLB - loopExp);
    EpsilonZ = Power(10, FirstDerivEpsilon - loopExp);
    Inter1.lo(Str) = Power(10, Inter1Epsilon - loopExp - 1);
    epsiSmax = Power(10, -5 - loopExp);

*    FlowsheetOpt_CEOS9b2_U.savePoint = 1;

    if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms

$include ../%ProbObjSpecInitFile%

      solve FlowsheetOpt_CEOS9b2_U using NLP minimizing Z9;
    else
*      solve FlowsheetOpt_CEOS9b2_R using NLP minimizing Z9;
    );
    T.l(Str) = Tscaled.l(Str)*Tref
    display T.l, ComplPen, epsiSmax, Inter0;
    ChckCmpl(0);

*$include ../ProcsFiles/RemoveUnusedEquipment.gms

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = GAMSCHK;
    );
    FlowsheetOpt_CEOS9b2_U.optfile = 0;
*    FlowsheetOpt_CEOS9b2_R.optfile = 0;

  );
);

if(ProblemTrayByTray eq 0,
  option NLP = convertd;

  if(HeatIntegModel eq 1,
    FlowsheetOpt_CEOS9b2_U.optfile = 1;
*    FlowsheetOpt_CEOS9b2_U.savepoint = 0;
    solve FlowsheetOpt_CEOS9b2_U using NLP minimizing Z9;
  else
*    FlowsheetOpt_CEOS9b2_R.optfile = 1;
*    FlowsheetOpt_CEOS9b2_R.savepoint = 0;
*    solve FlowsheetOpt_CEOS9b2_R using NLP minimizing Z9;
  );
);


  if(HeatIntegModel eq 1,
    solve FlowsheetOpt_CEOS9b2_U using NLP minimizing Z9;
  else
*    solve FlowsheetOpt_CEOS9b2_R using NLP minimizing Z9;
  );
  ChckCmpl(0);

  loop(Ed1,
    StgEd1.fx(Ed1) = round(StgEd1.l(Ed1));
  );

  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = GAMSCHK;
  );
*  FlowsheetOpt_CEOS9b2_U.savepoint = 1;
  FlowsheetOpt_CEOS9b2_U.optfile = 0;
  if(HeatIntegModel eq 1,
    solve FlowsheetOpt_CEOS9b2_U using NLP minimizing Z9;
  else
*    solve FlowsheetOpt_CEOS9b2_R using NLP minimizing Z9;
  );
  ChckCmpl(0);

  if(HeatIntegModel eq 1,
    Sec3Stat('modelstat') = FlowsheetOpt_CEOS9b2_U.modelstat;
    Sec3Stat('solvestat') = FlowsheetOpt_CEOS9b2_U.solvestat;
  else
*    CEOSStat('modelstat') = FlowsheetOpt_CEOS9b2_R.modelstat;
*    CEOSStat('solvestat') = FlowsheetOpt_CEOS9b2_R.solvestat;
  );

  execute_unload 'Sec3_CEOSFlowsheetResults.gdx';
  if(ProblemTrayByTray eq 0,
    option NLP = convertd;
*  FlowsheetOpt_CEOS9b2_R.optfile = 1;
*  FlowsheetOpt_CEOS9b2_R.savepoint = 0;
    FlowsheetOpt_CEOS9b2_U.optfile = 1;
*    FlowsheetOpt_CEOS9b2_U.savepoint = 0;

    if(HeatIntegModel eq 1,
      solve FlowsheetOpt_CEOS9b2_U using NLP minimizing Z9;
    else
*    solve FlowsheetOpt_CEOS9b2_R using NLP minimizing Z9;
    );
  );

*execute_unload 'matdata_connect.gdx', RealStr, InactiveStr, ShdwStr, NotShdwStr, InGnrlE, OutLGnrlE, OutVGnrlE, InOneGnrlE, FlashLiq, FlashVap, VapStr, LiqStr, HIZone, HIMap, CEqp, HEqp ;
);

if(ProblemCEOSFlowsheetOpt eq 2,
  execute_loadpoint 'Sec3_CEOSFlowsheetResults.gdx';
);

Sec3Stat('timer') = timeElapsed - timer;
timer = timeElapsed;

*4444444444444444444444444444444444*
************************************
********** BEGIN SECTION 4 *********
************************************
*4444444444444444444444444444444444*

SectionSwitch = 4;
$include ../%ProbSecSpecInitFile%

Variables
  LghtRec(Clmn)
  HvyRec(Clmn)
  Dist2Feed(Clmn);

*$include ../ProcsFiles/ClmnProcessing.gms


*EpsilonZ = FirstDerivEpsilon;
*EpsilonA = Inter1Epsilon;

*ComplPen = 1;

$include ../ProcsFiles/RemoveUnusedEquipment.gms

$include ../ProcsFiles/ConnectMESHwBypass3.gms
$include ../ProcsFiles/InitAndSolveMESHwBypass3.gms


Model FlowsheetOpt_TrayByTrayb_U /FlowsheetOpt_CEOS9b2_U - EqCsdCEOS - BubbleDewEqns - EqCascade - ObjSec3, TrayEqnsWithBypass, ObjSec4/;

* EqPRel1PReb
*Model FlowsheetOpt_TrayByTrayb_R /FlowsheetOpt_CEOS9b2_R - EqCsdCEOS - BubbleDewEqns - EqCascade - EqPRel1PReb, TrayEqns/;

Model FlowsheetOpt_TrayByTrayb2_U /FlowsheetOpt_CEOS_Base - EqCsdCEOS - BubbleDewEqns - EqCascade, HeatIntegUpdated, TrayEqnsWithBypass, ObjSec4/;
* EqPRel1PReb

*Model FlowsheetOpt_TrayByTrayb2_R /FlowsheetOpt_CEOS9b2_R - EqCsdCEOS - BubbleDewEqns - EqCascade  - EqPRel1PReb - obj9, TrayEqns, obj9b/;

if(ProblemTrayByTray eq 1,

$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%

  F.fx(Str)$(InactiveStr(Str) OR ManualRemove(Str)) = 0;
  Fc.fx(Str, Comp)$(InactiveStr(Str) OR ManualRemove(Str)) = 0;

  if(FirstSolverTray eq 1,
    option NLP = GAMSCHK;
  elseif FirstSolverTray eq 2,
    option NLP = CONOPT;
    FlowsheetOpt_TrayByTrayb_U.optfile = 1;
*    FlowsheetOpt_TrayByTrayb_R.optfile = 1;
  else
    option NLP = IPOPTH;
  );

*  FlowsheetOpt_TrayByTrayb_U.savepoint = 1;
*  FlowsheetOpt_TrayByTrayb_R.savepoint = 1;

eff.lo(Trays) = 0;

epsiSmax = 1E-3;

*option iterlim = 0;

$include ../InitFiles/ThrmEBubDew.gms
$include ../InitFiles/InitGnrlE.gms

saveAlpha = alpha;

*FlowsheetOpt_TrayByTrayb_U.savepoint = 2;

for(loopCount = 1 to 1 by 1,
  for(loopExp = 0 to 1 by 1,
* change back to 0 to 2 later
    ComplPen = Power(10, loopExp + StartCPPower);
    EpsilonZ = Power(10, FirstDerivEpsilon - loopExp);
    Inter1.lo(Str) = Power(10, Inter1Epsilon - loopExp);
    epsiSmax = Power(10, -4 - loopExp);

    if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
$include ../%ProbObjSpecInitFile%

*      FlowsheetOpt_TrayByTrayb_U.iterlim = 0;
*      solve FlowsheetOpt_TrayByTrayb_U using NLP minimizing Z9;
*      FlowsheetOpt_TrayByTrayb_U.iterlim = %GlobalIterlim%;
      solve FlowsheetOpt_TrayByTrayb_U using NLP minimizing Z9;
    else
*      solve FlowsheetOpt_TrayByTrayb_R using NLP minimizing Z9;
    );
    T.l(Str) = Tscaled.l(Str)*Tref;
    NumTrays.l(Ed1) = SUM(Trays$TrayCasc(Ed1,Trays), eff.l(Trays));
    display T.l, NumTrays.l, CascSize, ComplPen, epsiSmax;
    ChckCmpl(0);
$include ../ProcsFiles/RemoveUnusedEquipment.gms
$include ../ProcsFiles/AdjustMESHwBypass2.gms

    if(loopExp eq 1,
      alpha = -0.10;
    );

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = GAMSCHK;
    );
    FlowsheetOpt_TrayByTrayb_U.optfile = 0;
*    FlowsheetOpt_TrayByTrayb_R.optfile = 0;

  );
);

  solve FlowsheetOpt_TrayByTrayb_U using NLP minimizing Z9;

  if(HeatIntegModel eq 1,
    Sec4Stat('modelstat') = FlowsheetOpt_TrayByTrayb_U.modelstat;
    Sec4Stat('solvestat') = FlowsheetOpt_TrayByTrayb_U.solvestat;
  else
*    TrayByTrayStat('modelstat') = FlowsheetOpt_TrayByTrayb_R.modelstat;
*    TrayByTrayStat('solvestat') = FlowsheetOpt_TrayByTrayb_R.solvestat;
  );

  execute_unload 'Sec4_MESH_FlowsheetResults.gdx';

  if(ProblemHeatChecker eq 0,
    option NLP = convertd;
    FlowsheetOpt_TrayByTrayb_U.optfile = 1;
*    FlowsheetOpt_TrayByTrayb_U.savepoint = 0;
    solve FlowsheetOpt_TrayByTrayb_U using NLP minimizing Z9;
  );
);
* End if

*$include ../ProcsFiles/ClmnProcessing.gms
*$include ../ProcsFiles/FinalProcessing.gms

* May need to toggle on
*execute_unload %matout%;

*if(ProblemTrayByTray eq 1,
*execute_unload 'matdata_tray.gdx', F, Fc, Xc, T, Tscaled, P, ZEOS, H, S, Trays, TrayCasc, Ed1, InTrL, OutTrL, InTrV, OutTrV, Kt;

*execute_unload 'matdata_connect.gdx', RealStr, InactiveStr, TrayStr, ShdwStr, NotShdwStr, InGnrlE, OutLGnrlE, OutVGnrlE, InOneGnrlE, FlashLiq, FlashVap, VapStr, LiqStr, HIZone, HIMap, CEqp, HEqp ;

*);

if(ProblemTrayByTray eq 2,
  F.l(Str) = 0;
  Fc.l(Str,Comp) = 0;
  execute_loadpoint 'Sec4_MESH_FlowsheetResults.gdx';
);

$include ../ProcsFiles/RemoveUnusedEquipment.gms

Sec4Stat('timer') = timeElapsed - timer;
timer = timeElapsed;

*5555555555555555555555555555555555*
************************************
********** BEGIN SECTION 5 *********
************************************
*5555555555555555555555555555555555*


display "Beginning Section 5";
SectionSwitch = 5;
$include ../%ProbSecSpecInitFile%

$include ../HeatIntegrationChecker4.gms

if(ProblemHeatChecker eq 1,
  execute_unload 'Sec5_HeatCheckerResults.gdx';
elseif ProblemHeatChecker eq 2,
  execute_loadpoint 'Sec5_HeatCheckerResults.gdx';
);

Sec5Stat('timer') = timeElapsed - timer;
timer = timeElapsed;

*6666666666666666666666666666666666*
************************************
********** BEGIN SECTION 6 *********
************************************
*6666666666666666666666666666666666*

display "Beginning Section 6";

SectionSwitch = 6;

******************************************************************
***** Reoptimize with decomposed general heat exchange units *****
******************************************************************

* Restore RealStr and fEOSStr
RealStr(OldRealStr) = yes;
fEOSStr(OldfEOSStr) = yes;

* Remove splitter outlets from fEOSStr
loop(Sptr,
   fEOSStr(Str)$OutSptr(Sptr,Str) = no;
);

* Reinitialize CEOS advanced variables
$include ../InitFiles/Init_CEOS_Adv.gms

*display RealStr, fEOSStr, InactiveStr;

* Remove valves from InactiveGnrlE
InactiveGnrlE(Valve)$NotInactiveGnrlE(Valve) = no;
InactiveGnrlE(Flsh)$NotInactiveGnrlE(Flsh) = no;

* Add sub heat exchange units to heat integration formulation
loop(HIZone,
  loop(GnrlE$HIMap(GnrlE, HIZone),
    loop(SHtEx$GnrlMap(GnrlE,SHtEx),
      HIMap(SHtEx,HIZone) = yes;
    );
  );
);

* Remove decomposed units... may not be needed
loop(HIZone,
  HIMap(GnrlE, HIZone)$ActGnrl(GnrlE) = no;
);

display ActGnrl, HIMap, InactiveGnrlE, NotInactiveGnrlE;

$include ../ProcsFiles/HeatIntegrationProcessing.gms
CEqp(CoolingSHtEx) = yes;
HEqp(HeatingSHtEx) = yes;


* Reset bounds
$include ../InitFiles/Bounds.gms
$include ../%ProbSpecificBounds%
F.fx(Str)$InactiveStr(Str) = 0;
Fc.fx(Str,Comp)$InactiveStr(Str) = 0;

$include ../%ProbSecSpecInitFile%

*ZEOS.lo(Str) = 0;

* Need to generalize this
*F.lo('S12') = 0;

alpha = saveAlpha;

* Remove HRATPlusBD for HRAT

*aaEOS.up(Str) = 100;
*aaEOS.lo(Str) = -100;

* Ensure no vapor flow out of the condenser subunits

*loop(TCond,
*  loop(Str$OutLGnrlE(TCond,Str),
*    loop(SHtEx$OutLGnrlE(SHtEx,Str),
*      F.fx(Str2)$OutVGnrlE(SHtEx,Str2) = 0;
*    );
*  );
*);

* display OutLGnrlE, OutVGnrlE;

*option iterlim = 0;

$ontext
Equation
  EqHeatPen2;

EqHeatPen2..
  heatPen =e= SUM(SHtEx$GnrlMap('HX15',SHtEx), Qin(SHtEx))*SUM(SHtEx$GnrlMap('HX16',SHtEx), Qout(SHtEx));
$offtext

Model FlowsheetOpt_TrayByTrayb_Decomp_U /FlowsheetOpt_TrayByTrayb2_U - ObjSec4 - EqPRel1PReb, EqPresSHtEx, EqualDeltaT, ObjSec6/;

*option NLP = IPOPTH;
option NLP = GAMSCHK;

*option iterlim = 0;

*RelaxThermo(Str) = no;

*$include ../InitFiles/Init_CEOS_Basic4.gms

*display H.l;

* This is important... Otherwise H is initialized incorrectly
* $include ../InitFiles/Init_CEOS_Adv.gms

*display H.l;

$include ../InitFiles/ThrmEBubDew.gms

*$set matout3 "'matdata_heatchckr2.gdx', F, Tscaled, H, Qin, Qout, InHtExSub, OutLHtExSub, CoolingSHtEx, HeatingSHtEx, GnrlMap, QsZ, QwZ"

*display H.l;

if(ProblemTrayByTrayDecomp > 0,

display RealStr, fEOSStr, InactiveStr;

***** For debugging
FlowsheetOpt_TrayByTrayb_Decomp_U.savepoint = 2;

***** For debugging

for(loopCount = 1 to 1 by 1,
  for(loopExp = 0 to 2 by 1,
    ComplPen = Power(10, loopExp + StartCPPower);
    EpsilonZ = Power(10, FirstDerivEpsilon - loopExp);
    Inter1.lo(Str) = Power(10, Inter1Epsilon - loopExp);
    epsiSmax = Power(10, -4 - loopExp);

    if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
$include ../%ProbObjSpecInitFile%

      FlowsheetOpt_TrayByTrayb_Decomp_U.iterlim = 0;
      solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;
      FlowsheetOpt_TrayByTrayb_Decomp_U.iterlim = %GlobalIterlim%;

      solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;
    else
*      solve FlowsheetOpt_TrayByTrayb_R using NLP minimizing Z9;
    );
    T.l(Str) = Tscaled.l(Str)*Tref;
    NumTrays.l(Ed1) = SUM(Trays$TrayCasc(Ed1,Trays), eff.l(Trays));
    display T.l, NumTrays.l, CascSize, ComplPen, epsiSmax;
    ChckCmpl(0);

$include ../ProcsFiles/RemoveUnusedEquipment.gms
$include ../ProcsFiles/AdjustMESHwBypass2.gms

*    option iterlim = 0;

    if(loopExp eq 1,
      alpha = -0.10;
*      FlowsheetOpt_TrayByTrayb_Decomp_U.savepoint = 1;
    );

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = GAMSCHK;
    );
    FlowsheetOpt_TrayByTrayb_Decomp_U.optfile = 0;
*    FlowsheetOpt_TrayByTrayb_R.optfile = 0;

  );
);

*solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;

*FlowsheetOpt_TrayByTrayb_Decomp_U.optfile = 0;
*option NLP = CONOPT
*option iterlim = 0;
option NLP = GAMSCHK;
solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;

execute_unload 'Sec6_MESH_HtExDecomp_Results.gdx';

Sec6Stat('modelstat') = FlowsheetOpt_TrayByTrayb_Decomp_U.modelstat;
Sec6Stat('solvestat') = FlowsheetOpt_TrayByTrayb_Decomp_U.solvestat;

FlowsheetOpt_TrayByTrayb_Decomp_U.optfile = 1;
*FlowsheetOpt_TrayByTrayb_Decomp_U.savepoint = 0;
option NLP = CONVERTD;
solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;

*execute_unload %matout3%;

);

Sec6Stat('timer') = timeElapsed - timer;

execute_unload %matout%;

*option iterlim = 0;
FlowsheetOpt_TrayByTrayb_Decomp_U.optfile = 0;
*option NLP = IPOPTH;
*solve FlowsheetOpt_TrayByTrayb_Decomp_U using NLP minimizing Z9;


*if(ProblemTrayByTray eq 1,
*execute_unload 'matdata_tray_final.gdx', F, Fc, Xc, T, Tscaled, P, ZEOS, H, S, Trays, TrayCasc, Ed1, InTrL, OutTrL, InTrV, OutTrV, Kt;
*execute_unload 'matdata_connect_final.gdx', RealStr, InactiveStr, TrayStr, ShdwStr, NotShdwStr, InGnrlE, OutLGnrlE, OutVGnrlE, InOneGnrlE, FlashLiq, FlashVap, VapStr, LiqStr, HIZone, HIMap, CEqp, HEqp ;
*);

