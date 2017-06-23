***************************************
***** InitAndSolveMESHwBypass.gms *****
***************************************

* This file analyzes the Edmister cascades, essembles a MESH model with bypass,
* initializes the new model and solves it, all automatically.

**************************
***** Initialization *****
**************************

* In this section the cascade enterance/exit pressures, temperatures and flowrates
* are extracted from the Edmister model. Linear interpolation is then used to
* initialize these variables insider the cascade. The new streams internal to the
* cascade are then added to the sets that control thermodynamic equations.
* Finally the thermodynamic property calculations are initialized for the new
* streams.

***** Initialize Flows, Temperatures and Pressures *****

Parameters TBot, TTop, PBot, PTop, FTopV(Comp), FTopL(Comp), FBotV(Comp), FBotL(Comp);

T.l(Str) = Tscaled.l(Str)*Tref;

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
  count = 0;
*  count = 1;

  loop(Trays$TrayCasc(Ed1,Trays),
* Initialize "initially active" trays (no bypass - 100% efficiency)
    if(InitActive(Trays),
      Fc.l(Str,Comp)$OutTrV(Trays,Str) = count*(FTopV(Comp) - FBotV(Comp))/(InitCascSize(Ed1)) + FBotV(Comp);
*      Fc.l(Str,Comp)$InTrV(Trays,Str) = count*(FTopV(Comp) - FBotV(Comp))/(InitCascSize(Ed1) + 1) + FBotV(Comp);
      Fc.l(Str,Comp)$OutTrL(Trays,Str) = count*(FTopL(Comp) - FBotL(Comp))/(InitCascSize(Ed1)) + FBotL(Comp);

      T.l(Str)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = count*(TTop - TBot)/(InitCascSize(Ed1)) + TBot;

*      T.l(Str)$InTrV(Trays,Str) = count*(TTop - TBot)/(InitCascSize(Ed1) + 1) + TBot;
*      T.l(Str)$OutTrL(Trays,Str) = count*(TTop - TBot)/(InitCascSize(Ed1) + 1) + TBot;

      P.l(Str)$OutTrV(Trays,Str) = (PTop + PBot)/2;
      P.l(Str)$OutTrL(Trays,Str) = (PTop + PBot)/2;
      count = count + 1;

      eff.l(Trays) = 1;
    else
* Initialize "initially inactive" trays (complete bypass - 0% efficiency)
      Fc.l(Str,Comp)$OutTrV(Trays,Str) = FTopV(Comp);
      Fc.l(Str,Comp)$OutTrL(Trays,Str) = FTopL(Comp);

      T.l(Str)$OutTrV(Trays,Str) = TTop;
      T.l(Str)$OutTrL(Trays,Str) = TTop;

      P.l(Str)$OutTrV(Trays,Str) = PTop;
      P.l(Str)$OutTrL(Trays,Str) = PTop;

      eff.l(Trays) = 0;
    );

    F.l(Str)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = SUM(Comp, Fc.l(Str,Comp));
    Xc.l(Str,Comp)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = Fc.l(Str,Comp) / max(F.l(Str),1E-4);
  );

*  NumTrays.l(Ed1) = count - 1;
  NumTrays.l(Ed1) = count;

  loop(Trays$TrayCasc(Ed1,Trays),
      Xc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Xc.l(Str2,Comp));
      Xc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Xc.l(Str2,Comp));

      T.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), T.l(Str2));
      P.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), P.l(Str2));

      T.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), T.l(Str2));
      P.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), P.l(Str2));

*    if(InitActive(Trays),
      Fc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Fc.l(Str2,Comp));
      Fc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Fc.l(Str2,Comp));
*    else
*      Fc.l(Str,Comp)$(ETrV(Trays,Str) OR ETrL(Trays,Str)) = 0;
*    );
  );
);

display Fc.l, Xc.l, InitActive;

F.l(Str)$TrayStr(Str) = SUM(Comp, Fc.l(Str,Comp));

if(FixEffForInit eq 1,
  eff.fx(Trays) = eff.l(Trays);
);

Tscaled.l(Str) = T.l(Str)/Tref;

*display TBot, TTop, PBot, PTop, FTopV, FTopL, FBotV, FBotL;
display T.l, Tscaled.l, F.l;

***** Consider New Streams for Thermo Calcs *****

FlashLiq(Str)$(TrayStr(Str) AND (TrayLiqStr(Str) OR TrayLiqEquil(Str))) = yes;
FlashVap(Str)$(TrayStr(Str) AND (TrayVapStr(Str) OR TrayVapEquil(Str))) = yes;

fEOSStr(TrayStr) = yes;
RealStr(TrayStr) = yes;

HCalc(Str)$TrayStr(Str) = yes;
PhiCalc(Str)$TrayStr(Str) = yes;

***** Remove Bubble and Dew Point Shadow Streams *****
loop(Ed1,
  DewPoint(Str)$(OutVEd1(Ed1,Str)) = no;
  BubPoint(Str)$(OutLEd1(Ed1,Str)) = no;
);

$include ../ProcsFiles/ProcessShadowStreams.gms

*LiqStr(LiqStrAdvOnly) = no;
*VapStr(VapStrAdvOnly) = no;
*fEOSStr(AdvOnly) = no;
*RealStr(AdvOnly) = no;

***** Initialize Basic CEOS Thermo Calcs *****

$include ../InitFiles/Bounds.gms

*F.l(TrayStr) = SUM(Comp, Fc.l(TrayStr, Comp));
*Xc.l(Str,Comp)$TrayStr(Str) = Fc.l(Str,Comp)/max(F.l(Str), 1E-4);

aEOS.l(Str,j)$TrayStr(Str) = Power(1 + fw(j)*(1-sqrt(T.l(Str)/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
bmEOS.l(Str)$TrayStr(Str) = max( SUM(j, Xc.l(Str,j)*bEOS(j)), epsi) ;
amEOS.l(Str)$TrayStr(Str) = max( SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2))*(1-kEOS(j,j2)))), epsi);

*amEOS.l(Str) = SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)
*    *sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2)+epsi*(Power(aEOS.l(Str,j),2) + Power(aEOS.l(Str,j2),2)))*(1-kEOS(j,j2)))) ;

bbEOS.l(Str)$TrayStr(Str) = max( bmEOS.l(Str)*P.l(Str)/(8.314E-2*T.l(Str)), epsi) ;
aaEOS.l(Str)$TrayStr(Str) = amEOS.l(Str)*P.l(Str)/Power(8.314E-2*T.l(Str),2) ;

* TODO: Consider replacing with Init_CEOS_Basic.gms if this step frequently fails
ZEOS.l(Str)$(TrayStr(Str) AND (TrayLiqStr(Str) OR TrayLiqEquil(Str))) = bbEOS.l(Str) + 0.001;
ZEOS.l(Str)$(TrayStr(Str) AND (TrayVapStr(Str) OR TrayVapEquil(Str))) = 0.95;

V.l(Str)$TrayStr(Str) = ZEOS.l(Str)*R*T.l(Str)/max(P.l(Str),1E-6);

*T.fx(Str) = T.l(Str);
*Tscaled.lo(Str) = max(T.l(Str)*0.85, Tmin)/Tref;
*Tscaled.up(Str) = min(T.l(Str)*1.15, Tmax)/Tref;
*P.fx(Str) = P.l(Str);

*Xc.fx(Str,Comp) = Xc.l(Str,Comp);
*Xc.lo(Str,Comp) = 0;
*Xc.up(Str,Comp) = 1;

Xc.lo(Str,Comp)$TrayStr(Str) = min(0.95*Xc.l(Str,Comp),1);
Xc.up(Str,Comp)$TrayStr(Str) = max(1.05*Xc.l(Str,Comp),0);

phiEOS.up(TrayStr,Comp) = 50;
*phiEOS.up(Str,Comp) = Inf;

*Inter1.lo(Str) = 1E-5;
*EpsilonZ = 1E-5;

****************************************
***** Advanced CEOS Initialization *****
****************************************

* After initializing ZEOS with a good guess (something near 0.95 for vapor streams,
* something near 0.001 for liquid streams) a small optimization problem is solved
* to help further initialize CEOS variables. In this initialization problem
* stream temperature, pressure and composition are fixed.

* The following code formulates and solves the CEOS initialization problem.

InitStr1(Str) = yes;
InitStr2(Str) = no;
$include ../InitFiles/InitInitEqns.gms


*option iterlim = 0;

* Idea: Interpolation is not working...

if(ProblemTrayByTray eq 1,
  option NLP = GAMSCHK;
*   option NLP = CONOPT;
*  option NLP = IPOPTH;
  Solve InitCEOSTrays minimizing Z using NLP;
);

* With the basic CEOS variables initialized, the advanced CEOS equations
* (departure functions, etc) are initialized for the new tray streams.

Kt.l(Trays,Comp) = 1;
phiEOS.l(TrayStr,Comp) = 1;

*$include ../IncludeFiles/Init_CEOS_Adv.gms

HIG.l(Str) = 1E-3*SUM(J, Xc.l(Str,j)*(Tref**5 *CpIG('5',J)/5*(Tscaled.l(Str)**5 - 1) + Tref**4 *CpIG('4',J)/4*(Tscaled.l(Str)**4 - 1)
         + Tref**3 *CpIG('3',J)/3*(Tscaled.l(Str)**3 - 1) + Tref**2 *CpIG('2',J)/2*(Tscaled.l(Str)**2 - 1) + Tref*CpIG('1',J)*(Tscaled.l(Str) - 1)));

SIG.l(fEOSStr) = SUM(J, Xc.l(fEOSStr,J)*(Tref**4 *CpIG('5',J)/4*(Tscaled.l(fEOSStr)**4 - 1) + Tref**3 *CpIG('4',J)/3*(Tscaled.l(fEOSStr)**3 - 1)
         + Tref**2 *CpIG('3',J)/2*(Tscaled.l(fEOSStr)**2 - 1) + Tref*CpIG('2',J)*(Tscaled.l(fEOSStr) - 1)
         + CpIG('1',J)*log(Tscaled.l(fEOSStr)))) - Rsi*log(P.l(fEOSStr)/Pref);

***** Intermediate Variables *****
dadT.l(Str)$fEOSStr(Str) = -(R/2)*sqrt(omegaA/T.l(Str))
         *SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*(1-kEOS(j,j2))*
         (fw(j2)*sqrt(aEOS.l(Str,j)*Tc(j2)/Pc(j2)) + fw(j)*sqrt(aEOS.l(Str,j2)*Tc(j)/Pc(j)))));

delta.l(Str,j)$fEOSStr(Str) = 2*sqrt(aEOS.l(Str,j))*Sum(j2, Xc.l(Str,j2)*sqrt(aEOS.l(Str,j2))*(1 - kEOS(j, j2)))/amEOS.l(Str);

bRatio.l(Str,j)$fEOSStr(Str) = (Tc(j)/Pc(j))/max(Sum(j2, Xc.l(Str, j2)*Tc(j2)/Pc(j2)),epsi);

Inter1.l(Str)$fEOSStr(Str) = max(ZEOS.l(Str) - bbEOS.l(Str), Inter1.lo(Str));

Inter2.l(Str)$fEOSStr(Str) = (ZEOS.l(Str) - bbEOS.l(Str))/max(ZEOS.l(Str), 1E-8);

Inter3.l(Str)$fEOSStr(Str) = max((2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS + Inter0))/
         (2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS - Inter0)), Inter3.lo(Str));

***** Departure Functions *****

H.l(Str)$(fEOSStr(Str) AND HCalc(Str)) = HIG.l(Str) + 1E-1*(T.l(Str)*dadT.l(Str)
         - amEOS.l(Str))*log(Inter3.l(Str))/(bmEOS.l(Str)*Inter0)
         + 1E-3*Rsi*T.l(Str)*(ZEOS.l(Str)-1);

S.l(Str)$(fEOSStr(Str) AND SCalc(Str)) = SIG.l(Str)
         + 100*(R*log(Inter2.l(Str)) + R*log(ZEOS.l(Str)*Pref/P.l(Str))
         - 1/(bmEOS.l(Str)*Inter0)*dadT.l(Str)*log(Inter3.l(Str)));

phiEOS.l(Str,j)$(fEOSStr(Str) AND PhiCalc(Str)) = exp(bRatio.l(Str,j)*(ZEOS.l(Str)-1)-log(Inter1.l(Str))
         + aaEOS.l(Str)/(bbEOS.l(Str)*Inter0)*(bRatio.l(Str,j) - delta.l(Str,j))*(log(Inter3.l(Str)))) ;

display phiEOS.l;

*option iterlim = 0;

* Idea: Interpolation is not working...

if(ProblemTrayByTray eq 1,
*  option NLP = GAMSCHK;
   option NLP = CONOPT;
*  option NLP = IPOPTH;

  InitCEOSTrays2.savepoint = 1;

  Solve InitCEOSTrays2 minimizing Z using NLP;
);

* After initializing the basic CEOS variables, the stream pressures, temperatures
* and compositions must be unfixed.

Tscaled.lo(Str) = Tmin/Tref;
Tscaled.up(Str) = Tmax/Tref;
Tscaled.lo(FeedStr) = 300/Tref;

P.lo(Str) = Pmin;
P.up(Str) = Pmax;
Xc.lo(Str,Comp) = 0;
Xc.up(Str,Comp) = 1;

Kt.l(Trays,Comp) = 1;

display Kt.l, HIG.l, SIG.l, H.l, S.l;

display fEOSStr, ZEOS.l, LiqStr, VapStr, FlashLiq, FlashVap;

* After the advanced CEOS variables are initialized, another CEOS initialization
* problem is formulated and solved. In this new problem the difference between
* the flowrates from the Edmister cascade (from the previous full ASU solve) and
* the new tray-by-tray model is minimized. Inlet compositions are fixed
* for this final initialization problem.

Sets
  InOrOut(Str)           Streams in or out of the tray-by-tray cascades   / /;


loop(Ed1,
  InOrOut(Str)$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
  InOrOut(Str)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = yes;

* Fix inlet compositions
*  Xc.fx(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = Xc.l(Str,Comp);
);

*** Bounds *****
* F.lo(TrayStr) = 0.01;
F.lo(TrayStr) = 0;

options NLP = GAMSCHK;

betaT.l(Trays) = 1;

InitStr1(Str) = no;
InitStr1(InOrOut) = yes;
InitStr2(Str) = no;

$include ../InitFiles/InitInitEqns.gms

*option iterlim = 0;

eff.fx(Trays) = eff.l(Trays);

if(ProblemTrayByTray eq 1,
  option NLP = GAMSCHK;
*  solve ASU_Init_TrayByTray1 using NLP minimizing Z;
);

* option iterlim = 10;

eff.lo(Trays) = 0;
eff.up(Trays) = 1;

effTarget(Trays) = eff.l(Trays);

Kt.l(Trays,Comp) = SUM(Str2$ETrL(Trays,Str2), phiEOS.l(Str2,Comp))
         /max(SUM(Str$ETrV(Trays,Str), phiEOS.l(Str,Comp)),1E-5);

if(ProblemTrayByTray eq 1,
  option NLP = GAMSCHK;
*  option NLP = IPOPTH;
  ASU_Init_TrayByTray2.savepoint = 1;
  solve ASU_Init_TrayByTray2 using NLP minimizing Z;

*  option nlp = convert;
*  ASU_Init_TrayByTray2.optfile = 1;
*  solve ASU_Init_TrayByTray2 using NLP minimizing Z;

);

display fEOSStr, InactiveStr, InactiveGnrlE, RealStr;

T.l(Str) = Tscaled.l(Str)*Tref;

* Calculate difference between Edmister and initialized tray-by-tray model

Parameter DiffFc(Str,Comp);

loop(Ed1,
  DiffFc(Str,Comp)$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = Fc.l(Str,Comp)
         - FcTarget(Str,Comp);
  Xc.lo(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = 0;
  Xc.up(Str,Comp)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)) = 1;
);

display DiffFc;

* The tray-by-tray model is now initialized. Time for some optimization!

