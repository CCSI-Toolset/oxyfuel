**************************************
***** InitAndSolveTrayByTray.gms *****
**************************************

* This file analyzes the Edmister cascades, essembles a tray-by-tray model,
* initializes the new model and solves it, all automatically.

**************************
***** Initialization *****
**************************

* In this section the cascade enterance/exit pressures, temperatures and flowrates
* are extracted from the Edmister model. Linear interpolation is then used to
* initialize these variables also the cascade. The new streams internal to the
* cascade are then added to the sets that control thermodynamic equations.
* Finally the thermodynamic property calculations are initialized for the new
* streams.

***** Initialize Flows, Temperatures and Pressures *****

Parameters TBot, TTop, PBot, PTop, FTopV(Comp), FTopL(Comp), FBotV(Comp), FBotL(Comp);

loop(Ed1,
* Extract temperatures, pressures and flowrates
  loop(Str$OutLEd1(Ed1,Str),
    TBot = T.l(Str);
    PBot = P.l(Str);
    FBotL(Comp) = Fc.l(Str,Comp);
  );

  loop(Str$InVEd1(Ed1,Str),
    FBotV(Comp) = Fc.l(Str,Comp);
  );

  loop(Str$OutVEd1(Ed1,Str),
    TTop = T.l(Str);
    PTop = P.l(Str);
    FTopV(Comp) = Fc.l(Str,Comp);
  );

  loop(Str$InLEd1(Ed1,Str),
    FTopL(Comp) = Fc.l(Str,Comp);
  );

* Initialize the tray temperatures, pressures and flowrates using linear
* interpolation
  count = 1;
  loop(Trays$(TrayCasc(Ed1,Trays) AND NOT BotTrays(Ed1,Trays)),
    Fc.l(Str,Comp)$InTrV(Trays,Str) = count*(FTopV(Comp) - FBotV(Comp))/(CascSize(Ed1) + 1) + FBotV(Comp);
    Fc.l(Str,Comp)$OutTrL(Trays,Str) = count*(FTopL(Comp) - FBotL(Comp))/(CascSize(Ed1) + 1) + FBotL(Comp);
    T.l(Str)$(InTrV(Trays,Str) OR OutTrL(Trays,Str)) = count*(TTop - TBot)/(CascSize(Ed1) + 1) + TBot;
    P.l(Str)$(InTrV(Trays,Str) OR OutTrL(Trays,Str)) = (PTop + PBot)/2;
    count = count + 1;
  );
);

***** Consider New Streams for Thermo Calcs *****

LiqStr(Str)$(TrayStr(Str) AND TrayLiqStr(Str)) = yes;
VapStr(Str)$(TrayStr(Str) AND TrayVapStr(Str)) = yes;

fEOSStr(TrayStr) = yes;
RealStr(TrayStr) = yes;

HCalc(Str)$TrayStr(Str) = yes;
PhiCalc(Str)$TrayStr(Str) = yes;

***** Remove Bubble and Dew Point Shadow Streams *****
LiqStr(LiqStrAdvOnly) = no;
VapStr(VapStrAdvOnly) = no;
fEOSStr(AdvOnly) = no;

***** Initialize Basic CEOS Thermo Calcs *****

Tscaled.l(Str) = T.l(Str)/Tref;

F.l(TrayStr) = SUM(Comp, Fc.l(TrayStr, Comp));
Xc.l(Str,Comp)$TrayStr(Str) = Fc.l(Str,Comp)/max(F.l(Str), 1E-4);

aEOS.l(Str,j)$TrayStr(Str) = Power(1 + fw(j)*(1-sqrt(T.l(Str)/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
bmEOS.l(Str)$TrayStr(Str) = max( SUM(j, Xc.l(Str,j)*bEOS(j)), epsi) ;
amEOS.l(Str)$TrayStr(Str) = max( SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2))*(1-kEOS(j,j2)))), epsi);

*amEOS.l(Str) = SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)
*    *sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2)+epsi*(Power(aEOS.l(Str,j),2) + Power(aEOS.l(Str,j2),2)))*(1-kEOS(j,j2)))) ;

bbEOS.l(Str)$TrayStr(Str) = max( bmEOS.l(Str)*P.l(Str)/(8.314E-2*T.l(Str)), epsi) ;
aaEOS.l(Str)$TrayStr(Str) = amEOS.l(Str)*P.l(Str)/Power(8.314E-2*T.l(Str),2) ;

ZEOS.l(Str)$(TrayLiqStr(Str) AND TrayStr(Str)) = bbEOS.l(Str) + 0.001;
ZEOS.l(Str)$(TrayVapStr(Str) AND TrayStr(Str)) = 0.95;

V.l(Str)$TrayStr(Str) = ZEOS.l(Str)*R*T.l(Str)/max(P.l(Str),1E-6);

T.fx(Str) = T.l(Str);
Tscaled.lo(Str) = T.l(Str)*0.95/Tref;
Tscaled.up(Str) = T.l(Str)*1.05/Tref;
P.fx(Str) = P.l(Str);
Xc.fx(Str,Comp) = Xc.l(Str,Comp);

****************************************
***** Advanced CEOS Initialization *****
****************************************

* After initializing ZEOS with a good guess (something near 0.95 for vapor streams,
* something near 0.001 for liquid streams) a small optimization problem is solved
* to help further initialize CEOS variables. In this initialization problem
* stream temperature, pressure and composition are fixed.

* The following code formulates and solves the CEOS initialization problem.

Model InitCEOSTrays /InitObj2, BasicCEOSEqns/;

if(ProblemTrayByTray eq 1,
  option NLP = GAMSCHK;
  Solve InitCEOSTrays minimizing Z using NLP;
);

* After initializing the basic CEOS variables, the stream pressures, temperatures
* and compositions must be unfixed.

T.lo(Str) = Tmin;
T.up(Str) = Tmax;
Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;
P.lo(Str) = Pmin;
P.up(Str) = Pmax;
Xc.lo(Str,Comp) = 0;
Xc.up(Str,Comp) = 1;

* Strange stuff happens when this is an up not fx
T.fx(FeedStr) = 300;

* With the basic CEOS variables initialized, the advanced CEOS equations
* (departure functions, etc) are initialized for the new tray streams.

Kt.l(Trays,Comp) = 1;
phiEOS.l(TrayStr,Comp) = 1;

$include ./InitFiles/Init_CEOS_Adv.gms

display phiEOS.l;

Kt.l(Trays,Comp) = SUM(Str2$OutTrL(Trays,Str2), phiEOS.l(Str2,Comp))
         /max(SUM(Str$OutTrV(Trays,Str), phiEOS.l(Str,Comp)),1E-5);


*Kt.l(Trays,Comp) = 1;

display Kt.l, HIG.l, SIG.l, H.l, S.l;

display fEOSStr, ZEOS.l, LiqStr, VapStr;

* After the advanced CEOS variables are initialized, another CEOS initialization
* problem is formulated and solved. In this new problem the difference between
* the flowrates from the Edmister cascade (from the previous full ASU solve) and
* the new tray-by-tray model is minimized. Inlet compositions are fixed
* for this final initialization problem.

Parameters
  FcEdTarget(Str,Comp);

Sets
  InOrOut(Str)           Streams in or out of the tray-by-tray cascades   / /;


loop(Ed1,
  InOrOut(Str)$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
  InOrOut(Str)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = yes;

* Fix inlet compositions
  Xc.fx(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = Xc.l(Str,Comp);
);

* Extract flowrate targets
FcEdTarget(Str,Comp)$InOrOut(Str) = Fc.l(Str,Comp);


display FcEdTarget;

*** Bounds *****
F.lo(TrayStr) = 0.01;

Equation
  EqInitObj4;

EqInitObj4..
  Z =e= SUM(Comp, SUM(Ed1, SUM(Str$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str) OR InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)),
                 (Fc(Str,Comp) - FcEdTarget(Str,Comp))*
                 (Fc(Str,Comp) - FcEdTarget(Str,Comp)) )));

Z.l = SUM(Comp, SUM(Str$(InOrOut(Str)),
                 (Fc.l(Str,Comp) - FcEdTarget(Str,Comp))*
                 (Fc.l(Str,Comp) - FcEdTarget(Str,Comp)) ));

options NLP = GAMSCHK;

Model ASU_Init_TrayByTray1 /TrayEqns - EqTrayEnrgBal,CubicEOSEqns,EqStreams,EqInitObj4/;

Model ASU_Init_TrayByTray2 /TrayEqns,CubicEOSEqns,EqStreams - EqTscaled,EqInitObj4/;

Tscaled.l(Str) = T.l(Str)/Tref;
Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;

if(ProblemTrayByTray eq 1,
  solve ASU_Init_TrayByTray2 using NLP minimizing Z;
);

* Calculate difference between Edmister and initialized tray-by-tray model

Parameter DiffFc(Str,Comp);

loop(Ed1,
  DiffFc(Str,Comp)$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = Fc.l(Str,Comp)
         - FcEdTarget(Str,Comp);
  Xc.lo(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = 0;
  Xc.up(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = 1;
);

display DiffFc, FcEdTarget;

* The tray-by-tray model is now initialized. Time for some optimization!
