***************************************
***** HeatIntegrationChecker4.gms *****
***************************************

* This file takes existing heat integration designs and recalculates composite
* curves with several intermediate points for each heat exchanger. The file is
* organized into three main parts:

* 0. Load models, set declarations
* 1. Parsing existing heat integration design
* 2. Prepare optimization model
* 3. Solve thermodynamic questions posed as an optimization problem
* 4. Generate composite curves in post processing

* This new file addresses General Heat Exchange Units, not just heat exchangers

$onmacro
$onempty
$onexpand

****************************
********** Part 0 **********
****************************

*Scalar
*  epsi      small number for zero situation                         /1E-4/
*  SelectCEOSModel indicate which CEOS thermo model                  /3/;

*$include ../SpecFiles/StreamsComponents_complex.gms
*$include ../SpecFiles/FlowsheetTopology_26b_complex.gms
*$include ../ProcsFiles/FlowsheetTopologyProcessing.gms
*$include ../SpecFiles/Thermo_Data.gms
*$include ../ModelFiles/StreamModel.gms
*$include ../ModelFiles/ThrmEModel.gms
*$include ../ModelFiles/CondenserModel.gms
*$include ../ModelFiles/CommonThermoModel.gms
*$include ../ModelFiles/CEOSThermoModel.gms

*$include ../ModelFiles/SetupPostProcessing.gms

Sets
  InHtExSub(SHtEx, Str)            Inlet streams                   / /
  OutVHtExSub(SHtEx, Str)          Vapor outlet streams            / /
  OutLHtExSub(SHtEx, Str)          Liquid outlet streams           / /
  CoolingSHtEx(SHtEx)              Cooling subunits                / /
  HeatingSHtEx(SHtEx)              Heating subunits                / /
  SLiqStr(Str)                     Liquid streams                  /L1*L500/
  SVapStr(Str)                     Vapor streams                   /V1*V500/
  OutGnrlEStr(Str)                                                 / /;

Alias (SHtEx, SHtEx2);

display F.l;

****************************
********** Part 1 **********
****************************

if(1 eq 0,

* Initialize values
  Qin.l(GnrlE) = 0;
  Qout.l(GnrlE) = 0;
*  QTCond.l(TCond) = 0;
  F.l(Str) = 0;
  Fc.l(Str,Comp) = 0;
  T.l(Str) = 0;
  P.l(Str) = 0;
  ZEOS.l(Str) = 0;

* Bounds
T.lo(Str) = 65;
T.up(Str) = 400;
Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;

* Load existing heat integration solution
* execute_loadpoint "matdata_tray.gdx";
* execute_loadpoint "ASU_TrayByTrayb_U_p.gdx";

*display T.l, P.l, F.l;

);

*****
* Step A: Determine which heat exchangers, condensers and reboilers are active

Sets
  ActHtEx(HtEx)                                    / /
  ActReb(PReb)                                     / /
  ActCnd(TCond)                                    / /
  ActGnrl(GnrlE)
  NotInactiveGnrlE(GnrlE)                          / /
  ActSHtEx(SHtEx)                                  / /
  ActStr(Str)                                      / /
*  TwoPhase(HtEx)                                   / /
*  LiqOnly(HtEx)                                    / /
*  VapOnly(HtEx)                                    / /
  FixFc(Str)                                       / /;

  ActHtEx(HtEx)$(SUM(Str$InHtEx(HtEx,Str), F.l(Str)) > 0 OR Qin.l(HtEx) + Qout.l(HtEx) > 0) = yes;
  ActReb(PReb) = yes$(SUM(Str$InPReb(PReb,Str), F.l(Str)) > 0 OR Qin.l(PReb) + Qout.l(PReb) > 0);
  ActCnd(TCond) = yes$(SUM(Str$InTCond(TCond,Str), F.l(Str)) > 0 OR Qout.l(Tcond) > 0);

  ActGnrl(HtEx)$ActHtEx(HtEx) = yes;
  ActGnrl(PReb)$ActReb(PReb) = yes;
  ActGnrl(TCond)$ActCnd(TCond) = yes;

  display ActHtEx, ActReb, ActCnd, ActGnrl;

*****
* Step B: Allocate subunits to heat exchangers, condensers and reboilers
* Step C: Allocate streams to subunits
* Step D: Initialize subunits
* Step E: Fix stream properties for big unit entrance and exit

Scalar
  TotNumSubUnits                                 /0/
  SubUnitsPerUnit                                /4/
  n            Counter                           /0/
  nOffset      Offset from MESH model            /0/
  FLout                                          /0/
  FVout                                          /0/
  FLin                                           /0/
  FVin                                           /0/
  ZEOSLin                                        /0/
  ZEOSLout                                       /0/
  ZEOSVin                                        /0/
  ZEOSVout                                       /0/
  dZEOSL                                         /0/
  dZEOSV                                         /0/
  slackSmall                                     /0.1/
  slackLarge                                     /10/;

sL.up(Str) = 100;
sV.up(Str) = 100;

TotNumSubUnits = card(SHtEx);

Sets
  GnrlMap(GnrlE,SHtEx)                            / / ;

Parameters
*  Ttarget(Str)                                   / /
*  Ftarget(Str)                                   / /

  dT(GnrlE)                                      / /
  FcLin(Comp)                                / /
  FcVin(Comp)                                / /
  FcLout(Comp)                               / /
  FcVout(Comp)                               / /
  dFcV(Comp)                                 / /
  dFcL(Comp)                                 / /
  Fwght(Str)                                 / /;

* n = ord(AllocT(Trays)) + 1;

nOffset = 350;
n = nOffset;

*Preserve T for streams not involved in heat integration
Ttarget(Str)$(Not ActStr(Str)) = Tscaled.l(Str);

*display InGnrlE, OutVGnrlE, OutLGnrlE;

* Heat exchangers providing cooling
loop(GnrlE$ActGnrl(GnrlE),
* New Step - record F target
  Ftarget(Str)$(OutLGnrlE(GnrlE,Str) OR InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)) = F.l(Str);

  loop(Str$(OutLGnrlE(GnrlE,Str) OR InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)),
    if(F.l(Str) > 1E-4,
      Fwght(Str) = 100;
    else
      Fwght(Str) = 0;
    );
  );

* Step B
  GnrlMap(GnrlE,SHtEx)$(ord(SHtEx) > n AND ord(SHtEx) < n + SubUnitsPerUnit + 1) = yes;

* Step C
  loop(SHtEx$(GnrlMap(GnrlE,SHtEx) AND ord(SHtEx) < n + SubUnitsPerUnit),
    OutLHtExSub(SHtEx,SLiqStr)$(ord(SHtEx) eq ord(SLiqStr)) = yes;
    OutVHtExSub(SHtEx,SVapStr)$(ord(SHtEx) eq ord(SVapStr)) = yes;
  );

  loop(SHtEx2$(ord(SHtEx2) > n AND ord(SHtEx2) < n + SubUnitsPerUnit + 1),
* Want in_i+1 = out_i
    InHtExSub(SHtEx, Str)$(GnrlMap(GnrlE,SHtEx) AND ord(SHtEx) > n + 1 AND (ord(SHtEx) eq ord(SHtEx2) + 1) AND (OutVHtExSub(SHtEx2,Str) OR OutLHtExSub(SHtEx2,Str))) = yes;
  );

* Exceptions: inlet for first subunit, outlet for last subunit
  InHtExSub(SHtEx, Str)$( GnrlMap(GnrlE,SHtEx) AND InGnrlE(GnrlE, Str) AND ord(SHtEx) eq n + 1 ) = yes;
  OutVHtExSub(SHtEx, Str)$( GnrlMap(GnrlE,SHtEx) AND OutVGnrlE(GnrlE, Str) AND ord(SHtEx) eq n + SubUnitsPerUnit ) = yes;
  OutLHtExSub(SHtEx, Str)$( GnrlMap(GnrlE,SHtEx) AND OutLGnrlE(GnrlE, Str) AND ord(SHtEx) eq n + SubUnitsPerUnit ) = yes;

* Step D
  dT(GnrlE) = ( SUM(Str$OutLGnrlE(GnrlE, Str), Tscaled.l(Str)) - SUM(Str2$InGnrlE(GnrlE, Str2), Tscaled.l(Str2))/SUM(Str2$InGnrlE(GnrlE,Str2), 1) )/(SubUnitsPerUnit);

  FLout = SUM(Str$OutLGnrlE(GnrlE,Str), F.l(Str)) + 0;
  FVout = SUM(Str$OutVGnrlE(GnrlE,Str), F.l(Str)) + 0;

  FLin = SUM(Str$(InGnrlE(GnrlE,Str) AND (LiqStr(Str) OR FlashLiq(Str))), F.l(Str)) + 0;
  FVin = SUM(Str$(InGnrlE(GnrlE,Str) AND (VapStr(Str) OR FlashVap(Str))), F.l(Str)) + 0;

  FcLout(Comp) = SUM(Str$OutLGnrlE(GnrlE,Str), Fc.l(Str,Comp)) + 0;
  FcVout(Comp) = SUM(Str$OutVGnrlE(GnrlE,Str), Fc.l(Str,Comp)) + 0;

  FcLin(Comp) = SUM(Str$(InGnrlE(GnrlE,Str) AND (LiqStr(Str) OR FlashLiq(Str))), Fc.l(Str,Comp)) + 0;
  FcVin(Comp) = SUM(Str$(InGnrlE(GnrlE,Str) AND (VapStr(Str) OR FlashVap(Str))), Fc.l(Str,Comp)) + 0;

  dFcL(Comp) = (FcLout(Comp) - FcLin(Comp))/SubUnitsPerUnit;
  dFcV(Comp) = (FcVout(Comp) - FcVin(Comp))/SubUnitsPerUnit;

* Is this used? - now it is
  ZEOSLin = SUM(Str$(InGnrlE(GnrlE,Str) AND (LiqStr(Str) OR FlashLiq(Str))), ZEOS.l(Str)) + 0;
  ZEOSVin = SUM(Str$(InGnrlE(GnrlE,Str) AND (VapStr(Str) OR FlashVap(Str))), ZEOS.l(Str)) + 0;

  ZEOSLout = SUM(Str$OutLGnrlE(GnrlE,Str), ZEOS.l(Str)) + 0;
  ZEOSVout = SUM(Str$OutVGnrlE(GnrlE,Str), ZEOS.l(Str)) + 0;

  dZEOSL = (ZEOSLout - ZEOSLin)/SubUnitsPerUnit;
  dZEOSV = (ZEOSVout - ZEOSVin)/SubUnitsPerUnit;

  loop(SHtEx$GnrlMap(GnrlE,SHtEx),
    Ttarget(Str)$InHtExSub(SHtEx, Str) = SUM(Str2$InGnrlE(GnrlE, Str2), Tscaled.l(Str2))/SUM(Str2$InGnrlE(GnrlE,Str2), 1) + dT(GnrlE)*(ord(SHtEx) - n - 1);
    Tscaled.l(Str)$InHtExSub(SHtEx, Str) = Ttarget(Str);

    Fc.l(Str,Comp)$OutLHtExSub(SHtEx, Str) = FcLin(Comp) + dFcL(Comp)*(ord(SHtEx) - n);
    Fc.l(Str,Comp)$OutVHtExSub(SHtEx, Str) = FcVin(Comp) + dFcV(Comp)*(ord(SHtEx) - n);

    loop(Str$(OutVHtExSub(SHtEx, Str) OR OutLHtExSub(SHtEx, Str)),

      StrNow(Str2) = no;
      StrNow(Str) = yes;

      P.l(Str) = SUM(Str2$InHtExSub(SHtEx, Str2), P.l(Str2))/SUM(Str2$InHtExSub(SHtEx,Str2), 1);
      F.l(Str) = SUM(Comp, Fc.l(Str,Comp));

      if(F.l(Str) < epsi,
          Xc.l(Str,Comp) = (FcLin(Comp) + FcVin(Comp))/(FLin + FVin);
      else
        Xc.l(Str,Comp) = Fc.l(Str,Comp)/F.l(Str);
      );
    );

$ontext
    if( FLout > 0 AND FVout > 0,
* Two phase outlet
      ZEOS.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = 0.01;
      sL.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackSmall;
      ZEOS.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = 0.95;
      sV.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackSmall;
      Fc.l(Str,Comp)$(OutVHtExSub(SHtEx, Str) OR OutLHtExSub(SHtEx, Str)) = SUM(Str2$InHtExSub(SHtEx, Str2), Fc.l(Str2,Comp))/2;
*      TwoPhase(HtEx) = yes;

    elseif (FLin > 0 AND FVout > 0) OR (FVin > 0 and FLout > 0),
*     elseif 1 < 0,
* Two phase inlet, either liquid only or vapor only outlet

      if(FVout > 0,
        ZEOS.l(Str)$OutLHtExSub(SHtEx, Str) = ZEOSLin + dZEOSL*(ord(SHtEx) - n - 1);
        ZEOS.l(Str)$OutVHtExSub(SHtEx, Str) = 0.95;
        sV.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackSmall;
        sL.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackLarge;
      else
        ZEOS.l(Str)$OutVHtExSub(SHtEx, Str) = ZEOSVin + dZEOSV*(ord(SHtEx) - n - 1);
        ZEOS.l(Str)$OutLHtExSub(SHtEx, Str) = 0.01;
        sV.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackLarge;
        sL.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackSmall;
      );

      Fc.l(Str,Comp)$OutLHtExSub(SHtEx, Str) = FcLin(Comp) + dFcL(Comp)*(ord(SHtEx) - n - 1);
      Fc.l(Str,Comp)$OutVHtExSub(SHtEx, Str) = FcVin(Comp) + dFcV(Comp)*(ord(SHtEx) - n - 1);

    elseif FLout > 0,
* Liquid only outlet
*      ZEOS.l(Str)$((SLiqStr(Str) OR SVapStr(Str)) AND (OutLHtExSub(SHtEx, Str) OR OutVHtExSub(SHtEx, Str))) = 0.01;
      ZEOS.l(Str)$( (SLiqStr(Str) OR SVapStr(Str)) AND (OutLHtExSub(SHtEx, Str) OR OutVHtExSub(SHtEx, Str)))
*                 = ZEOSin + dZEOS*(ord(SHtEx) - n - 1);
*                 = ZEOSout;
                 = 0.01;
      sL.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackSmall;
      sL.up(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackSmall;
      sV.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackLarge;
*      LiqOnly(HtEx) = yes;
      Fc.l(Str,Comp)$OutVHtExSub(SHtEx, Str) = 0.01*SUM(Str2$InHtExSub(SHtEx, Str2), Fc.l(Str2,Comp));
      Fc.l(Str,Comp)$OutLHtExSub(SHtEx, Str) = 0.99*SUM(Str2$InHtExSub(SHtEx, Str2), Fc.l(Str2,Comp));

    else
* Vapor only outlet
*      ZEOS.l(Str)$((SLiqStr(Str) OR SVapStr(Str)) AND (OutLHtExSub(SHtEx, Str) OR OutVHtExSub(SHtEx, Str))) = 0.95;
      ZEOS.l(Str)$((SLiqStr(Str) OR SVapStr(Str)) AND (OutLHtExSub(SHtEx, Str) OR OutVHtExSub(SHtEx, Str)))
*                 = ZEOSin + dZEOS*(ord(SHtEx) - n - 1);
*                 = ZEOSin;
                  = 0.95;
*      VapOnly(HtEx) = yes;

      sL.l(Str)$(SLiqStr(Str) AND OutLHtExSub(SHtEx, Str)) = slackLarge;
      sV.l(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackSmall;
      sV.up(Str)$(SVapStr(Str) AND OutVHtExSub(SHtEx, Str)) = slackSmall;
      Fc.l(Str,Comp)$OutLHtExSub(SHtEx, Str) = 0.01*SUM(Str2$InHtExSub(SHtEx, Str2), Fc.l(Str2,Comp));
      Fc.l(Str,Comp)$OutVHtExSub(SHtEx, Str) = 0.99*SUM(Str2$InHtExSub(SHtEx, Str2), Fc.l(Str2,Comp));

    );
$offtext

*display Xc.l;


  );

* Step E
  P.fx(Str)$InGnrlE(GnrlE,Str) = P.l(Str);
*  Fc.fx(Str,Comp)$(InHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str) OR OutLHtEx(HtEx,Str)) = Fc.l(Str,Comp);
*  Fc.fx(Str,Comp)$InHtEx(HtEx,Str) = Fc.l(Str,Comp);
  FixFc(Str)$InGnrlE(GnrlE,Str) = yes;

* Fixing Tscaled seems to cause some problems
*  Tscaled.fx(Str)$(InHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str) OR OutLHtEx(HtEx,Str)) = T.l(Str)/Tref;

  Tscaled.lo(Str)$(InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = Tscaled.l(Str) - 2/Tref;
  Tscaled.up(Str)$(InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = Tscaled.l(Str) + 2/Tref;
*  Tscaled.fx(Str)$(InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = T.l(Str)/Tref;

  Ttarget(Str)$(InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = Tscaled.l(Str);

* Step F
  if( CEqp(GnrlE),
    CoolingSHtEx(SHtEx)$GnrlMap(GnrlE, SHtEx) = yes;
  elseif HEqp(GnrlE),
    HeatingSHtEx(SHtEx)$GnrlMap(GnrlE, SHtEx) = yes;
  );

* Finish loop
  n = n + SubUnitsPerUnit;
);

loop(GnrlE$ActGnrl(GnrlE),
  FixFc(Str)$(OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = no;
);

* Problem (F.l > 0) not introduced yet

*****
* Step G: "Deactivate" large units

NotInactiveGnrlE(ThrmE)$(NOT SHtEx(ThrmE) AND NOT InactiveGnrlE(ThrmE)) = yes;
NotInactiveGnrlE(TCond)$(NOT InactiveGnrlE(TCond)) = yes;
InactiveGnrlE(ThrmE)$(NOT SHtEx(ThrmE)) = yes;
InactiveGnrlE(TCond) = yes;


Qout.up(SHtEx)$CoolingSHtEx(SHtEx) = Inf;
Qin.up(SHtEx)$HeatingSHtEx(SHtEx) = Inf;

*****
* Step H: Save fEOSStr and reset

Sets
  OldRealStr(Str)  / /
  OldfEOSStr(Str)  / /;

OldRealStr(RealStr) = yes;
OldfEOSStr(fEOSStr) = yes;

RealStr(Str) = no;
fEOSStr(Str) = no;

display GnrlMap, ActStr, InHtExSub, OutLHtExSub, OutVHtExSub, Ttarget, dT;

*****
* Step I: Populate GnrlE connectivity sets using small units

* Thermo Equip Connectivity
loop(SHtEx,
  InGnrlE(SHtEx,Str)$(InHtExSub(SHtEx,Str)) = yes;
  OutVGnrlE(SHtEx,Str)$(OutVHtExSub(SHtEx,Str)) = yes;
  OutLGnrlE(SHtEx,Str)$(OutLHtExSub(SHtEx,Str)) = yes;

* Set stream pairs
  StrPairs(Str,Str2)$(OutVGnrlE(SHtEx,Str) AND OutLGnrlE(SHtEx,Str2)) = yes;
);

loop(SHtEx,
  loop(Str$(InGnrlE(SHtEx, Str) AND NOT InactiveStr(Str)),
    if(SUM(Str2$InOneGnrlE(SHtEx, Str2), 1) < 1,
      InOneGnrlE(SHtEx, Str) = yes;
    );
  );
);

* Initially populate InThrmEPres with all of the inlet streams for all thermo
* equipment.
InThrmEPres(ThrmE,Str)$(InGnrlE(ThrmE,Str)) = yes;

* Remove the liquid stream in each pair (leaving ThrmE A) if the
* vapor and liquid stream are both feeds for a different unit (feeds for ThrmE B)
loop(ThrmE,
  loop(Str$InGnrlE(ThrmE,Str),
    loop(Str2$InGnrlE(ThrmE,Str2),
      if(StrPairs(Str,Str2),
        InThrmEPres(ThrmE,Str2) = no;
      );
    );
  );
);

* Remove exit casacde stream from InThrmEPres(Flsh,Str)
* This prevents a redundent pressure equation
loop(Ed1,
  loop(Flsh,
    InThrmEPres(Flsh,Str)$(OutLEd1(Ed1,Str)) = no;
  );
);

loop(ThrmE,
  loop(Str$(OutVGnrlE(ThrmE,Str) OR OutLGnrlE(ThrmE,Str)),
*   Collect all of the outlet thrme streams
    OutGnrlEStr(Str) = yes;
  );
);

* Usage to remove recycles. Consider StrPairs(1,2). If 1 is an outlet for unit A,
* don't consider 2 as an inlet.

loop(ThrmE,
  loop(Str$(OutLGnrlE(ThrmE,Str) OR OutVGnrlE(ThrmE,Str)),
    loop(Str2$StrPairs(Str,Str2),
      InThrmEPres(ThrmE,Str2) = no;
    );
  );
);


* Populate SHtEx
loop(SHtEx$ActSHtEx(SHtEx),
  ActStr(Str)$(InHtExSub(SHtEx,Str) OR OutVHtExSub(SHtEx,Str) OR OutLHtExSub(SHtEx,Str)) = yes;
);

* Add sub heat exchangers back in
ActSHtEx(SHtEx)$(ord(SHtEx) < n + 1) = yes;
CalcGnrlE(SHtEx)$ActSHtEx(SHtEx) = yes;

Qin.fx(CoolingSHtEx) = 0;
Qout.fx(HeatingSHtEx) = 0;


****************************************
***** Model Declaration and Bounds *****
****************************************

****
* Step A: Declare models

Positive Variables
  liqPen
  vapPen

Variables
  Tdiff(Str)
  Fdiff(Str)
  slack(SHtEx)
  deltaT(GnrlE);

Variables
  Z;

Equations
  EqObj1
  EqObj2
  EqObj3
  EqObj4
  EqTdiff(Str)
  EqFdiff(GnrlE,Str)
  EqLiqPen2               Calculate liquid penalty - used for complementarities
  EqVapPen2               Calculate vapor penalty - used for complementarities
  EqPresSHtEx(SHtEx,Str,Str2)            No pressure drop in subunit HX
  EqOrderCoolingSHtEx    Ensure vapor stream decreases in cooling subunits
  EqOrderHeatingSHtEx    Ensure liquid stream decreases in heating subunits
  EqCalcTotalQin(GnrlE)
  EqCalcTotalQout(GnrlE)
  EqEqualQoutSHtEx(GnrlE,SHtEx)
  EqEqualQinSHtEx(GnrlE,SHtEx)
  EqCalcDeltaT(GnrlE,Str,Str2)
  EqEqualDeltaT(GnrlE,SHtEx,Str,Str2)
  EqCoolSHtEx(SHtEx,Str,Str2)
  EqHeatSHtEx(SHtEx,Str,Str2);

EqObj1..
  Z =e= SUM(Str$ActStr(Str), 100*Tdiff(Str)*Tdiff(Str) + sL(Str) + sV(Str) + 0.1*Fdiff(Str)*Fdiff(Str));

EqObj2..
  Z =e= SUM(Str$ActStr(Str), 100*Tdiff(Str)*Tdiff(Str) + 0.1*Fdiff(Str)*Fdiff(Str)) + 100*liqPen + 100*vapPen;

EqObj3..
  Z =e= SUM(Str$ActStr(Str), Fdiff(Str)*Fdiff(Str)) + 100*liqPen + 100*vapPen + 10*SUM(Str$FixFc(Str), SUM(Comp, Power(XcTarget(Str,Comp) - Xc(Str,Comp),2) ));

EqObj4..
*  Z =e= 1000*liqPen + 1000*vapPen + 0.01*SUM(ThrmE$(NOT InactiveGnrlE(ThrmE)), Power(1 - beta(ThrmE), 2));
Z =e= 1000*liqPen + 1000*vapPen;

EqTdiff(Str)$ActStr(Str)..
  Tdiff(Str) =e= Tscaled(Str) - Ttarget(Str);

EqFdiff(GnrlE,Str)$(ActGnrl(GnrlE) AND (OutLGnrlE(GnrlE,Str) OR InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)))..
  Fdiff(Str) =e= F(Str) - Ftarget(Str);

EqVapPen2..
 vapPen =e= sum(Str$(ActStr(Str) AND FlashVap(Str)), sV(Str)*F(Str));

EqLiqPen2..
 liqPen =e= sum(Str$(ActStr(Str) AND FlashLiq(Str)), sL(Str)*F(Str));

EqPresSHtEx(SHtEx,Str,Str2)$(InThrmEPres(SHtEx,Str) AND OutLHtExSub(SHtEx,Str2) AND ActSHtEx(SHtEx) AND NOT InactiveStr(Str))..
  P(Str) =e= P(Str2);

EqOrderCoolingSHtEx(SHtEx,Str,Str2)$(InHtExSub(SHtEx,Str) AND SVapStr(Str) AND OutVHtExSub(SHtEx,Str2))..
  F(Str) =g= F(Str2);

EqOrderHeatingSHtEx(SHtEx,Str,Str2)$(InHtExSub(SHtEx,Str) AND SLiqStr(Str) AND OutLHtExSub(SHtEx,Str2))..
  F(Str) =g= F(Str2);

EqCalcTotalQin(GnrlE)$(HEqp(GnrlE))..
  Qin(GnrlE) =e= SUM(SHtEx$GnrlMap(GnrlE, SHtEx), Qin(SHtEx));

EqCalcTotalQout(GnrlE)$(CEqp(GnrlE))..
  Qout(GnrlE) =e= SUM(SHtEx$GnrlMap(GnrlE, SHtEx), Qout(SHtEx));

EqEqualQinSHtEx(GnrlE,SHtEx)$(HeatingSHtEx(SHtEx) AND GnrlMap(GnrlE, SHtEx))..
  Qin(SHtEx) =e= Qin(GnrlE)/SubUnitsPerUnit + slack(SHtEx);

EqEqualQoutSHtEx(GnrlE,SHtEx)$(CoolingSHtEx(SHtEx) AND GnrlMap(GnrlE, SHtEx))..
  Qout(SHtEx) =e= Qout(GnrlE)/SubUnitsPerUnit + slack(SHtEx);

EqCalcDeltaT(GnrlE,Str,Str2)$((HEqp(GnrlE) OR CEqp(GnrlE)) AND ActGnrl(GnrlE) AND InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2))..
  deltaT(GnrlE) =e= Tscaled(Str) - Tscaled(Str2);

EqEqualDeltaT(GnrlE,SHtEx,Str,Str2)$(GnrlMap(GnrlE,SHtEx) AND InOneGnrlE(SHtEx,Str) AND OutLGnrlE(SHtEx,Str2))..
  Tscaled(Str) - Tscaled(Str2) =e= deltaT(GnrlE)/SubUnitsPerUnit;

EqCoolSHtEx(CoolingSHtEx,Str,Str2)$(InThrmEPres(CoolingSHtEx,Str) AND OutLGnrlE(CoolingSHtEx,Str2) AND NOT InactiveGnrlE(CoolingSHtEx))..
  Tscaled(Str) =g= Tscaled(Str2);

EqHeatSHtEx(HeatingSHtEx,Str,Str2)$(InThrmEPres(HeatingSHtEx,Str) AND OutLGnrlE(HeatingSHtEx,Str2) AND NOT InactiveGnrlE(HeatingSHtEx))..
  Tscaled(Str) =l= Tscaled(Str2);

Model EqOrderSHtEx /EqOrderCoolingSHtEx, EqOrderHeatingSHtEx/;

Model EqualQ /EqCalcTotalQin, EqCalcTotalQout, EqEqualQinSHtEx, EqEqualQoutSHtEx/;

*Model EqualDeltaT /EqCalcDeltaT, EqEqualDeltaT, EqCoolSHtEx, EqHeatSHtEx/;
*Model EqualDeltaT /EqEqualDeltaT, EqCoolSHtEx, EqHeatSHtEx/;

Model EqualDeltaT /EqEqualDeltaT/;

*Model HeatCheckerInit /EqStreams - EqTScaled, EqThrmE - EqEnrgBalGnrlE - EqVLEThrmE, BasicCEOSEqns, EqObj1, EqTdiff, EqLiqPen, EqVapPen, EqPresSHtEx/
Model HeatCheckerInit /EqStreams - EqTScaled, EqThrmE - EqEnrgBalGnrlE - EqVLEThrmE, BasicCEOSEqns, EqObj1, EqTdiff, EqFdiff, EqLiqPen2, EqVapPen2, EqPresSHtEx, EqCoolSHtEx, EqHeatSHtEx/
* EqOrderSHtEx

Model HeatCheckerInit2 /EqStreams - EqTScaled, EqThrmE - EqEnrgBalGnrlE - EqVLEThrmE, CubicEOSEqns, EqObj1, EqTdiff, EqFdiff, EqLiqPen2, EqVapPen2, EqPresSHtEx, EqCalcTotalQin, EqCalcTotalQout, EqualDeltaT/
* EqOrderSHtEx

Model InitSubunits /EqStreams - EqTScaled, EqThrmE, CubicEOSEqns, EqObj4, EqLiqPen2, EqVapPen2, EqPresSHtEx/;

Model HeatChecker /EqStreams - EqTScaled, EqThrmE, CubicEOSEqns, EqObj3, EqFdiff, EqLiqPen2, EqVapPen2, EqPresSHtEx, EqualDeltaT, EqCalcTotalQin, EqCalcTotalQout/;
* EqOrderSHtEx

* Display the automatically populated sets for debugging purposes.
display InGnrlE, OutVGnrlE, OutLGnrlE, StrPairs, InThrmEPres, OutGnrlEStr, InactiveGnrlE;

*****
* Step B: Establish bounds

EpsilonZ = 1E-8;
EpsilonA = 1E-8;

phiEOS.lo(fEOSStr,j) = epsi;
phiEOS.up(fEOSStr,j) = 1000;

* Set bounds... potential for divide by zero
bmEOS.lo(fEOSStr) = epsi;
bbEOS.lo(fEOSStr) = epsi;
amEOS.lo(fEOSStr) = epsi;

* Potential sqrt of a negative number
aEOS.lo(fEOSStr,Comp) = 0;

* Avoid sqrt of a negative number
aEOS.lo(fEOSStr,Comp) = epsi;

ZEOS.lo(Str) = epsi*epsi;
ZEOS.up(Str) = 5;

******************************************
***** Initialize subunit by subunit *****
******************************************


Scalars n2, n3;

Set FcFix2(Str), SkipZEOSInit(Str);

* Don't initialize ZEOS variables for inlet/outlet streams because they already have values
* The ZEOS initialization is sensitivity to rounding errors I suspect happen during GDX imports
SkipZEOSInit(Str) = no;
loop(ActGnrl,
  SkipZEOSInit(Str)$(InGnrlE(ActGnrl,Str) OR OutVGnrlE(ActGnrl,Str) OR OutLGnrlE(ActGnrl,Str)) = yes;
);

display Tscaled.l;

* n2 is the current "iteration"
* Trick: Don't need to initialize the exist streams of large units, as they are already in VLE
for(n2 = 1 to SubUnitsPerUnit - 1 by 1,

*****
* Step A: Reset by emptying sets

  ActStr(Str) = no;
  fEOSStr(Str) = no;
  RealStr(Str) = no;
  InactiveGnrlE(SHtEx) = yes;
  ActSHtEx(SHtEx) = no;

*****
* Step B: Determine active streams and subunits

* Find subunits corresponding to the current "iteration"
  loop(GnrlE,
    n3 = 1;
    loop(SHtEx$GnrlMap(GnrlE,SHtEx),
      if(n2 eq n3,
        InactiveGnrlE(SHtEx) = no;
        ActStr(Str)$(InGnrlE(SHtEx,Str) OR OutVGnrlE(SHtEx,Str) OR OutLGnrlE(SHtEx,Str)) = yes;
        FcFix2(Str)$InGnrlE(SHtEx,Str) = yes;
        ActSHtEx(SHtEx) = yes;
      );
      n3 = n3 + 1;
    );
  );

  display ActStr, ActSHtEx;

*****
* Step C: Activate thermo models for active streams

  fEOSStr(Str)$ActStr(Str) = yes;
  RealStr(Str)$ActStr(Str) = yes;
  HCalc(Str)$ActStr(Str) = yes;
  PhiCalc(Str)$ActStr(Str) = yes;
  FlashLiq(Str)$(ActStr(Str) AND SLiqStr(Str)) = yes;
  FlashVap(Str)$(ActStr(Str) AND SVapStr(Str)) = yes;

*****
* Step D: Initialize CEOS model

*$batinclude ../InitFiles/Init_CEOS_Basic4.gms "(ActStr(Str) AND NOT SkipZEOSInit(Str))"
$batinclude ../InitFiles/Init_CEOS_Basic4.gms "ActStr(Str)"
  display ZEOS.l;

$batinclude ../InitFiles/Init_CEOS_Adv2.gms "ActStr(Str)"

*****
* Step E: Fix inlet stream flowrates and pressures for subunits being initialized

* Added this rounding code on May 22nd
  loop(ActStr$(F.l(ActStr) < 1E-6),
    F.l(ActStr) = 0;
    Fc.l(ActStr,Comp) = 0;
  );


  Fc.fx(Str,Comp)$FcFix2(Str) = Fc.l(Str,Comp);
  P.fx(Str)$FcFix2(Str) = P.l(Str);

*****
* Step F: Fix temperature (instead of using temperature sparing constraint)
  Tscaled.fx(ActStr) = Tscaled.l(ActStr);

*****
* Step G: Initialize Qin and Qout

loop(GnrlE,
  Qin.l(GnrlE) = SUM(Str$(OutLGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)), F.l(Str)*H.l(Str)) - SUM(Str$InGnrlE(GnrlE,Str), F.l(Str)*H.l(Str));
  Qout.l(GnrlE) = SUM(Str$InGnrlE(GnrlE,Str), F.l(Str)*H.l(Str)) - SUM(Str$(OutLGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)), F.l(Str)*H.l(Str));
  deltaT.l(GnrlE) = SUM(Str$InOneGnrlE(GnrlE,Str), Tscaled.l(Str)) - SUM(Str$OutLGnrlE(GnrlE,Str), Tscaled.l(Str));
);

*****
* Step H: Solve NLP to initialize subunits

  if(ProblemHeatChecker > 0,

    option NLP = GAMSCHK;
    InitSubunits.holdfixed = 1;

*  InitSubunits.iterlim = 0;
*  Solve InitSubunits using NLP minimizing Z;
*  ChckCmpl(0);

    InitSubunits.iterlim = %GlobalIterlim%;
    Solve InitSubunits using NLP minimizing Z;
    ChckCmpl(0);

    if(InitSubunits.ModelStat eq 5,
      option NLP = MINOS;
      Solve InitSubunits using NLP minimizing Z;
      ChckCmpl(0);
    );
$ontext
  if(InitSubunits.ModelStat ne 2,
    option NLP = SNOPT;
    Solve InitSubunits using NLP minimizing Z;
    ChckCmpl(0);
  );

  if(InitSubunits.ModelStat ne 2,
    option NLP = IPOPTH;
    Solve InitSubunits using NLP minimizing Z;
    ChckCmpl(0);
  );
$offtext

    display Tscaled.l;
  );

);


***********************************
***** Finish Processing/Setup *****
*****          and            *****
*****   Solve HeatChecker     *****
***********************************

*****
* Step A: Reset ActSHtEx and assign all subunit streams to ActStr

ActSHtEx(SHtEx)$(ord(SHtEx) < n + 1) = yes;

loop(SHtEx$ActSHtEx(SHtEx),
  ActStr(Str)$(InHtExSub(SHtEx,Str) OR OutVHtExSub(SHtEx,Str) OR OutLHtExSub(SHtEx,Str)) = yes;
);

*****
* Step B: Activate thermo models for active streams

  fEOSStr(Str)$ActStr(Str) = yes;
  RealStr(Str)$ActStr(Str) = yes;

*****
* Step C: Remove all active subunits from InactiveGnrlE

InactiveGnrlE(SHtEx)$CoolingSHtEx(SHtEx) = no;
InactiveGnrlE(SHtEx)$HeatingSHtEx(SHtEx) = no;

*****
* Step D: Adjust fixed flowrates and temperatures

Fc.lo(Str,Comp)$FcFix2(Str) = 0;
Fc.up(Str,Comp)$FcFix2(Str) = 1000;

Fc.fx(Str,Comp)$FixFc(Str) = Fc.l(Str,Comp);

XcTarget(Str,Comp)$FixFc(Str) = Xc.l(Str,Comp);

*****
* Step E: Ensure Qin and Qout are initialized properly for the last set of subunits

$batinclude ../InitFiles/Init_CEOS_Adv2.gms "ActStr(Str)"

loop(GnrlE,
  Qin.l(GnrlE) = SUM(Str$(OutLGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)), F.l(Str)*H.l(Str)) - SUM(Str$InGnrlE(GnrlE,Str), F.l(Str)*H.l(Str));
  Qout.l(GnrlE) = SUM(Str$InGnrlE(GnrlE,Str), F.l(Str)*H.l(Str)) - SUM(Str$(OutLGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str)), F.l(Str)*H.l(Str));
  deltaT.l(GnrlE) = SUM(Str$InOneGnrlE(GnrlE,Str), Tscaled.l(Str)) - SUM(Str$OutLGnrlE(GnrlE,Str), Tscaled.l(Str));
);


*****
* Step F: Solve HeatChecker

if(ProblemHeatChecker > 0,
  HeatChecker.holdfixed = 1;
  HeatChecker.savepoint = 2;
*  HeatChecker.iterlim = 0;
*  Solve HeatChecker minimizing Z using NLP;

  HeatChecker.iterlim = %GlobalIterlim%;
  Solve HeatChecker minimizing Z using NLP;

  if(HeatChecker.ModelStat eq 5 OR Z.l > 1E-1,
    Solve HeatChecker minimizing Z using NLP;
  );

  if(HeatChecker.ModelStat eq 5 OR Z.l > 1E-1,
    option NLP = IPOPTH;
    Solve HeatChecker minimizing Z using NLP;
  );

  ChckCmpl(0);

  Sec5Stat('modelstat') = HeatChecker.modelstat;
  Sec5Stat('solvestat') = HeatChecker.solvestat;

);


$set matout2 "'matdata_heatchckr.gdx', F, Tscaled, H, Qin, Qout, InHtExSub, OutLHtExSub, CoolingSHtEx, HeatingSHtEx, GnrlMap"

if(ProblemHeatChecker > 0,
  execute_unload %matout2%;

  option NLP = CONVERTD;
  HeatChecker.optfile = 1;
  Solve HeatChecker minimizing Z using NLP;

  display ActSHtEx, InThrmEPres;
);

Z9.l = Z.l;

***************************************************************
***** Prepare results for integration with full flowsheet *****
***************************************************************

Tscaled.lo(Str) = Tmin/Tref;
Tscaled.up(Str) = Tmax/Tref;




