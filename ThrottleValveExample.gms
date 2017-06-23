* TOP OF FILE

* This is the "main" file for ASU synthesis and optimization. Include files are
* extensively used.
*
* The optimization problem is solved in four stages.
*
* 1. Simple thermodynamic ASU optimziation
*    In the first stage a simplified ASU model featuring correlation based
*    thermodynamic equations is optimized and heat integrated. This determines
*    reasonable temperatures, pressures, flows, etc for all of the ASU streams.                                                      )
*
* 2. Cubic CEOS model initialization
*    In the second stage stream properties are fixed, and the CEOS model
*    are initialized. This is done by solving three optimization problems.
*    During this stage the shadow streams used for bubble/dew point calculations
*    are "activated".
*
* 3. CEOS ASU model optimization with Edmister cascades
*    In the third stage the ASU with CEOS thermo is optimized and heat
*    integrated. The group method Edmister-based cascade approximation is used
*    in this model, which results in the number of theoretical trays being a
*    continous optimization variable. The optimization parameter is originally
*    solved with a large smoothing epsilon and loose group method bounds.
*    Epsilon is decrease and the bounds are tightened in subsequent solves. In
*    the final solves the number of stages in the cascades is rounded and fixed.
*
* 4. CEOS ASU model optimization with tray-by-tray cascades
*    To facilitate model validation, the Edmister cascade model is replaced with
*    a tray-by-tray model in the fourth stage. For these final solved the number
*    of trays remains fixed. Conversion from the group method to cascade model
*    is automated with dynamic sets in GAMS. Initialization is also automated
*    with loops.



* global options
$eolcom //
*option reslim = 6000;
option optcr = 1E-10; ;
option solprint = ON ;
option reslim = 3600;
$onmacro
$onempty
$onexpand

// end of global options

***** Stage Mode Settings *****
* These parameters control the "operational mode" for each stage.
Parameter
ThermoSwitch       Toggle between thermo methods                           /1/;
* Mode 1 - simple thermo (vapor pressure correlations)
* Mode 2 - CEOS thermo
* Should start as 1... don't change here


Scalars
Tmax      Max temperature (K) for the flowsheet                   /320/

Pmax      Max pressure (bar) for flowsheet                        /40/
Pmin      Min pressure (bar) for flowsheet                        /1.01/

* for zero flow in cascade
epsi      small number for zero situation                         /1E-4/

* Not important for this example
LBStgEd1  Lower bound on stages in cascade                       /2.0/
UBStgEd1  Upper bound on stages in cascade                       /50.0/

MaxDT     Maximium temperature difference for heat integration    /30/;
;

*------------------------**
* External Initialization
*------------------------**

Parameters
  NumStageInit   Number of initial stages in cascades
  PhiLBCEOS      Lower bound for phi in cascade model for CEOS problems
  UBASFact       Upper bound on absorption or stripping factor
  LBASFact       Lower bound on absorption or stripping factor
  Tmin           Min temperature (K) for the flowsheet
  Prune          Manually remove equipment before ASU_Simple?
  alpha          Small number used for heat integration phase changes
  HRATSimple     HRAT for ASU_Simple
  HRATCEOS       HRAT for ASU_CEOS and TrayByTray                        /2/
  O2RecLowDelta  Difference in lower bound for O2 recovery compared to CEOS (ASU with Simple thermo)
  O2RecLowCEOS   Lower bound for O2 recovery (ASU with CEOS thermo)
  ASUCEOSCompl       Value for complementarity penality for first ASU_CEOS solve
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
  LoadInitPoint Load initial flows-temps-pressures from file          /0/
  SimpleCscInitProb                                                   /0/
  CEOSCscInitProb                                                     /0/
  PruneConfig
  o2pureSpec
;

  NumStageInit = 10;
  UBASFact = 7;
  LBASFact = 0.001;
  Tmin = 60;
  Prune = 1;
  alpha = -0.15;
  HRATSimple = 6;
  O2RecLowDelta = 0;
  O2RecLowCEOS = 0.55;
  ASUCEOSCompl = 1;

  FirstSolverSimple = 1;
  FirstSolverCEOS = 1;
  FirstSolverTray = 1;

  HeatIntegModel = 1;

  PhiLBCEOS = -9;
  CascInterLB = -7;
  FirstDerivEpsilon = -7;
  Inter1Epsilon = -7;
  PlusStages = 15;
  FixEffForInit = 0;
  SelectCEOSModel = 3;
  PruneConfig = 4;

$include ../SpecFiles/StreamsComponents_example.gms

Set FeedStr(Str) /S1/
    ColumnFeedStr(Str) / /;

Parameters
  FeedFlow(Comp) /N2 = 0.78, O2 = 0.21, Ar = 0.01/
  ComplPen /100/
  epsiSmax /1E-4/
  loopExp /1/
  loopCount /1/;

$macro smmax(x) ( 0.5*sqrt(Power(x,2) + epsiSmax) + 0.5*(x) )


Scalars
  FeedT          /100/
  FeedP          /5/;

$include ../SpecFiles/Thermo_Data.gms

* Flowsheet Topology
$include ../SpecFiles/FlowsheetTopology_example.gms
$include ../ProcsFiles/FlowsheetTopologyProcessing.gms

HRAT('Z1') = HRATSimple;

************************************
***** Include flowsheet models *****
************************************

* Stream model
$include ../ModelFiles/StreamModel.gms

Parameters FeedTTmp /100/;

Fc.l(Str,Comp) = FeedFlow(Comp);
F.l(Str) = SUM(Comp, Fc.l(Str,Comp));
Xc.l(Str,Comp) = Fc.l(Str,Comp)/Max(F.l(Str),epsi);
Tscaled.l(Str) = FeedTTmp/Tref;
P.l(Str) = FeedP;

*********************************************
***** Continue include flowsheet models *****
*********************************************

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
$include ../ModelFiles/CEOSDewBubModel.gms

* Equipment models dependent on thermo
$include ../ModelFiles/CondenserModel.gms
$include ../ModelFiles/CascadeModel.gms
$include ../ModelFiles/EntropyBalModel.gms
$include ../ModelFiles/SplitterModel.gms

* Heat integration model
$include ../ModelFiles/HeatIntegrationModelCommon.gms
$include ../ModelFiles/HeatIntegrationModelUpdated.gms
*$include ../ModelFiles/HeatIntegrationModelRevised.gms

* Advanced distillation models - not needed for this example
* $include ../ModelFiles/TrayByTrayModel.gms
* $include ../ModelFiles/MESHwBypassModel3.gms

EpsilonZ = FirstDerivEpsilon;
EpsilonA = Inter1Epsilon;
Inter1.lo(Str) = 10**EpsilonA;

*--------------------------------*
*     Dummy variables
*--------------------------------*
Variables
Z ;

* Make sure energy balance for valve is active
CalcGnrlE(Valve) = yes;

$include ../ModelFiles/ObjFncsExample.gms
$include ../ModelFiles/SetupPostProcessing.gms


VapStr('S1') = yes;

$include ../InitFiles/Bounds.gms
$include ../InitFiles/BoundsExampleOnly.gms


*------------------------**
* Initialize Flows
*------------------------**

***** Option 1: Load From File *****
* In this option temperatures, pressures, compositions and flowrates
* are loaded from a file. This is analogous to starting from a known solution,
* such as results obtained from a commerical flowsheet simulator.

* Do nothing. This was already done above.

* The the throttle example, this option has been removed.

***** Option 2: Load a general guess ****
* Using this option the flowsheet is initialized to contain non-zero flows with
* some guess temperatures and pressures. In the current state, option 1 tends to
* be more effective. Refinement of the guesses in option 2 should improve
* performance.

Parameter FeedXc(Comp);

FeedXc(Comp) = FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));

Equation EqInitObj;

EqInitObj..
  Z =e= Sum(Str, (F(Str) - 1)*(F(Str) - 1)
         + 10*SUM(Comp, (Xc(Str,Comp) - FeedXc(Comp))*(Xc(Str,Comp) - FeedXc(Comp))));


Model MassBalEqns /EqTotMolStr, EqSplitFrac, EqTotCompMolBalEd1, EqCompMolBalGnrlE,
                  EqMolBalSptr, EqMolFracSptr/;

Model MassBal /MassBalEqns, EqInitObj/;

  MassBal.savepoint = 1;
  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = CONOPT;
  );

* EITHER: 1. Solve mass balance problem to initialize most streams with non-zero flows
* OR: 2. Don't overwrite initial point (F, Fc, Tscaled, P) loaded from file
  if(LoadInitPoint eq 0,

    Fc.l(Str,Comp) = FeedFlow(Comp);
    F.l(Str) = SUM(Comp, Fc.l(Str,Comp));
    Xc.l(Str,Comp) = Fc.l(Str,Comp)/Max(F.l(Str),epsi);

    Solve MassBal using NLP minimizing Z;

    T.l(Str) = FeedTTmp;
    P.l(Str) = FeedP;
    T.l(FeedStr) = 300;

    Tscaled.l(Str) = T.l(Str)/Tref;

  );

*------------------------**
* Prune Flowsheet
*------------------------**

Set
  ManualRemove(Str)      This set contains streams that should be manually removed. This forces certain flowsheet configurations        /  /;

* For the throttle example, no streams need to be pruned.

F.fx(ManualRemove) = 0;
Fc.fx(ManualRemove,Comp) = 0;

F.l(Str)$(NOT ManualRemove(Str)) = max(F.l(Str), 1E-3);

$include ../ProcsFiles/ManualRemoveEquipment.gms

display ManualRemove;

*********************************
*** Initialization and Bounds ***
*********************************

***** General Specs *****

T.lo(FeedStr) = 300;
T.up(FeedStr) = 350;
T.up(ColumnFeedStr) = 200;

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

***** Initialize Q in GnrlE *****

$include ../InitFiles/InitGnrlE.gms

***** Initialize Objective *****

Equation EqValveFixed(Valve,Str,Str2);

EqValveFixed(Valve,Str,Str2)$(InThrmEPres(Valve,Str) and OutLValve(Valve,Str2) AND NOT InactiveGnrlE(Valve))..
  P(Str) =g= P(Str2);

Model ThrottleExample /EqStreams,EqThrmE,EqFlsh,EqPReb,EqValveFixed,EqTCond,EqSpltr,EqObj/;

Model Throttle_Simple /ThrottleExample, SimpleThermo, EqSpltrSimpleThermo/;

Model Throttle_Simple_U /ThrottleExample - EqTscaled, SimpleThermo, EqSpltrSimpleThermo, EqHtEx, HeatIntegUpdated/;

**************************************************
***** Special Cascade Initialization Problem *****
**************************************************

Parameter
  FTarget(Str)
  FcTarget(Str,Comp)
  PTarget(Str)
  Ttarget(Str);

Variables
  FPen
  FcPen
  TPen
  PPen;

Equation
  EqInitCscObj2
  EqFPen
  EqFcPen
  EqTPen
  EqPPen;

EqFPen..
  FPen =e= SUM(Str$(NOT CscStr(Str)), (F(Str) - FTarget(Str))*(F(Str) - FTarget(Str)));

EqFcPen..
  FcPen =e= SUM(Comp, SUM(Str$CscStr(Str), (Fc(Str,Comp) - FcTarget(Str,Comp))*(Fc(Str,Comp) - FcTarget(Str,Comp))));

EqTPen..
  TPen =e= SUM(Str, (Ttarget(Str) - Tscaled(Str))*(Ttarget(Str) - Tscaled(Str)) );

EqPPen..
  PPen =e= SUM(Str, (PTarget(Str) - P(Str))*(PTarget(Str) - P(Str)));

EqInitCscObj2..
  Z =e= SUM(Str$CscStr(Str), sL(Str) + sV(Str)) + FPen
        + 0.1*FcPen
        + 0.1*TPen
        + 0.1*PPen;

FcTarget(Str,Comp) = Fc.l(Str,Comp);

Model InitCscObjs /EqFPen, EqFcPen, EqTPen, EqPPen, EqInitCscObj2/;

Model InitCscNewObj /MassBalEqns - EqTotCompMolBalEd1, EqSumFrac, SimpleThermo, EqCascade, EqCsdSimple, InitCscObjs/;
* EqSumFrac

Model InitMassBal /MassBalEqns, EqSumFrac, EqFcPen/;

* count = 0;

ComplPen = 10;

***** Initialize Heat Integration *****

if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
else
*$include ../InitFiles/Init_HtInteg_Revised.gms
);

sL.fx(LiqStr) = 0;
sV.fx(VapStr) = 0;

***** Initialize Slacks *****
  liqPen.l = sum(ThrmE, sum(Str$OutLGnrlE(ThrmE,Str), sL.l(Str)*F.l(Str)));
  vapPen.l = sum(ThrmE, sum(Str$OutVGnrlE(ThrmE,Str), sV.l(Str)*F.l(Str)));

  ComplPen = 10;
  epsiSmax = 1E-4;

$include ../InitFiles/ThrmEBubDew.gms
$include ../InitFiles/InitObjHelperEqnsExample.gms

  ComplPen = 1;

  if(FirstSolverSimple eq 1,
    option NLP = CONOPT;
  elseif FirstSolverSimple eq 2,
    option NLP = CONOPT;
*    Throttle_Simple_U_R.optfile = 1;
    Throttle_Simple_U.optfile = 1;
  elseif FirstSolverSimple eq 3,
    option NLP = IPOPTH;
  elseif FirstSolverSimple eq 4,
    option NLP = MINOS;
  elseif FirstSolverSimple eq 5,
    option NLP = SNOPT;
  else
    option NLP = CONOPT;
  );

*  Throttle_Simple_U_R.savepoint=1;
  Throttle_Simple_U.savepoint=1;

epsiSmax = 1E-3;

  Tscaled.l(Str) = T.l(Str)/Tref;
  Tscaled.lo(Str) = T.lo(Str)/Tref;
  Tscaled.up(Str) = T.up(Str)/Tref;

*option iterlim = 0;

for(loopCount = 1 to 1 by 1,
  for(loopExp = -1 to 0 by 1,

    ComplPen = Power(10, loopExp+2);

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
      solve Throttle_Simple_U using NLP minimizing Z;
    else
*     solve Throttle_Simple_U_R using NLP minimizing Z;
    );

    T.l(Str) = Tscaled.l(Str)*Tref;
    display T.l, ComplPen, epsiSmax;

    ChckCmpl(0);
*$include ../ProcsFiles/RemoveUnusedEquipment.gms

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = CONOPT;
    );
*    Throttle_Simple_U_R.optfile = 0;
    Throttle_Simple_U.optfile = 0;
  );
);

*$include ../ProcsFiles/AddUnusedEquipment.gms

  F.fx(ManualRemove) = 0;

* The following line dramatically impacts the solution
  F.l(Str)$(NOT ManualRemove(Str)) = max(F.l(Str), 10E-5);

* Tried this with initPoint4. No change when uncommented
* $include ../ProcsFiles/RemoveUnusedEquipment.gms

*******************************
********** Section 2 **********
*******************************
* CEOS Initialization

* Switch from Simple Thermo to CEOS Thermo Model
ThermoSwitch = 2;

Parameter XcFeed(Comp);

XcFeed(Comp) = FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));

Equation InitObj2;
InitObj2..
  Z =e= SUM(Str, sL(Str)) + SUM(Str, sV(Str));

Model InitCEOS /InitObj2, BasicCEOSEqns/;

$include ../InitFiles/PrepareInitCEOS2.gms

***** Initialize basic CEOS variables *****

$include ../InitFiles/Init_CEOS_Basic.gms


Z.l = sum(Str, sL.l(Str) + sV.l(Str));

RealStr(AdvOnly) = yes;
VapStr(VapStrAdvOnly) = yes;
LiqStr(LiqStrAdvOnly) = yes;

ZEOS.lo(Str) = 0;

Model InitCEOS1 /InitObj2, BasicCEOSEqns, EqSumFrac/;

*  option NLP = IPOPT;
  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = CONOPT;
  );
  Solve InitCEOS1 using NLP minimizing Z;

Model InitCEOS2 /InitObj2, CubicEOSEqns, EqSumFrac/;

$include ../InitFiles/Init_CEOS_Adv.gms

  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = CONOPT;
  );
*  option NLP = IPOPTH;
  Solve InitCEOS2 using NLP minimizing Z;


Parameter XcTarget(Str,Comp), Ttarget(Str);

XcTarget(Str,Comp) = Xc.l(Str,Comp);

Ttarget(Str) = Tscaled.l(Str);

Variables
  XcPen;

Equations InitObj3, EqXcPen;

EqXcPen..
  XcPen =e= SUM(Comp, SUM(Str2, SUM(Str$(BPStr(Str,Str2) OR DPStr(Str,Str2)),
           (XcTarget(Str,Comp) - Xc(Str,Comp))*(XcTarget(Str,Comp) - Xc(Str,Comp)) )));

InitObj3..
  Z =e= SUM(Str, sL(Str)) + SUM(Str, sV(Str))
         + 1000*XcPen + 0.1*SUM(Str$CscStr(Str), (Ttarget(Str) - Tscaled(Str))*(Ttarget(Str) - Tscaled(Str) ));

Model InitCEOS3 /InitObj3, EqXcPen, CubicEOSEqns, BubbleDewEqns, EqSumFrac, EqSpltrCEOSThermo, EqSpltr - EqMolBalSptr/;

Xc.lo(CscStr,Comp) = 0;
Xc.up(CscStr,Comp) = 1;

T.lo(CscStr) = T.l(CscStr) - 10;
T.up(CscStr) = T.up(CscStr) + 10;

Tscaled.lo(CscStr) = T.lo(CscStr)/Tref;
Tscaled.up(CscStr) = T.up(CscStr)/Tref;

  if(ExternalInit eq 1,
    option NLP = CONOPT;
  else
    option NLP = CONOPT;
  );
*  option NLP = IPOPTH;
  InitCEOS3.savepoint = 1;
  Solve InitCEOS3 using NLP minimizing Z;

*********************************
*** Initialization and Bounds ***
*********************************

***** General *****

$include ../InitFiles/Bounds.gms
$include ../InitFiles/BoundsExampleOnly.gms
epsi = 1e-6;

ComplPen = 100;

Xc.lo(Str,Comp) = 0;

T.lo(CscStr) = Tmin;

T.up(LiqStr) = smax(Comp, Tc(Comp)) + 5;

F.fx(LiqStrAdvOnly) = 1;
F.fx(VapStrAdvOnly) = 1;

Fc.l(ShdwStr,Comp) = Xc.l(ShdwStr,Comp);

Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;

F.fx(ManualRemove) = 0;

* Init HIG
$include ../InitFiles/Init_CEOS_Adv.gms

***** Cascade *****

***** Heat Exchangers *****

* Adjust zone approach temperatures
HRAT('Z1') = HRATCEOS;

***** Valves *****

* This line is neccissary with the JT specific valve model, but has no effect
* with the H(T,P) simple thermo correlation
CalcGnrlE(Valve) = yes;

Model Throttle_CEOS_Base /ThrottleExample, CubicEOSEqns, EqCsdCEOS, BubbleDewEqns, EqHtEx, EqSpltrCEOSThermo/;

*Model Throttle_CEOS_R /Throttle_CEOS_Base, HeatIntegRevised, obj9, EqCmprPwr/;
Model Throttle_CEOS_U /Throttle_CEOS_Base - EqTscaled, HeatIntegUpdated/;

$include ../InitFiles/InitGnrlE.gms

if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
else
*$include ../IncludeFiles/Init_HtInteg_Revised.gms
);

$include ../InitFiles/InitObjHelperEqnsExample.gms

  ComplPen = ASUCEOSCompl;

  PhiAEd1.lo(Ed1,Comp) = Power(10,PhiLBCEOS);
  PhiSEd1.lo(Ed1,Comp) = Power(10,PhiLBCEOS);

$include ../InitFiles/ThrmEBubDew.gms

if(Prune eq 1,
$include ../ProcsFiles/RemoveUnusedEquipment.gms
);

*  Throttle_CEOS_R.savepoint = 1;
  Throttle_CEOS_U.savepoint = 1;

  if(FirstSolverCEOS eq 1,
    if(ExternalInit eq 1,
      option NLP = CONOPT;
      Throttle_CEOS_U.optfile = 0;
    else
      option NLP = CONOPT;
    );
  elseif FirstSolverCEOS eq 2,
    option NLP = CONOPT;
    Throttle_CEOS_U.optfile = 1;
*    Throttle_CEOS_R.optfile = 1;
  else
    option NLP = IPOPTH;
  );

epsiSmax = 1E-3;

  Tscaled.l(Str) = T.l(Str)/Tref;
  Tscaled.lo(Str) = T.lo(Str)/Tref;
  Tscaled.up(Str) = T.up(Str)/Tref;

for(loopCount = 1 to 1 by 1,
*  for(loopExp = -1 to 1 by 1,
  for(loopExp = -1 to -1 by 1,
    ComplPen = Power(10, loopExp+1);
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

*    Throttle_CEOS_U.savePoint = 1;

    if(HeatIntegModel eq 1,
$include ../InitFiles/Init_HtInteg_Updated.gms
      Z.l = ComplPen*vapPen.l + ComplPen*liqPen.l
                 + Qw.l + Qs.l + SUM(FeedStr, CmprPwr.l(FeedStr));
      solve Throttle_CEOS_U using NLP minimizing Z;
    else
*      solve Throttle_CEOS_R using NLP minimizing Z;
    );
    T.l(Str) = Tscaled.l(Str)*Tref
    display T.l, ComplPen, epsiSmax, Inter0;
    ChckCmpl(0);

*$include ../ProcsFiles/RemoveUnusedEquipment.gms

    if(ExternalInit eq 1,
      option NLP = CONOPT;
    else
      option NLP = CONOPT;
    );
    Throttle_CEOS_U.optfile = 0;
*    Throttle_CEOS_R.optfile = 0;

  );

  if(HeatIntegModel eq 1,
    solve Throttle_CEOS_U using NLP minimizing Z;
  else
*    solve Throttle_CEOS_R using NLP minimizing Z;
  );
  ChckCmpl(0);
);

Variables
  LghtRec(Clmn)
  HvyRec(Clmn)
  Dist2Feed(Clmn);

*$include ../ProcsFiles/ClmnProcessing.gms
