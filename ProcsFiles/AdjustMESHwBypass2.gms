eff.l(Trays) = round(eff.l(Trays));
NumTrays.l(Ed1) = SUM(Trays$TrayCasc(Ed1,Trays), eff.l(Trays));

InitCascSize(Ed1) = round(NumTrays.l(Ed1));
CascSize(Ed1) = min(InitCascSize(Ed1) + PlusStages, MaxStagesPerCsc);

PrevAct(Ed1,Trays,Ordr) = no;

eff.up(Trays) = 1;

TopTraysOld(Ed1,Trays) = no;
TopTraysOld(Ed1,Trays)$TopTrays(Ed1,Trays) = yes;

loop(Ed1$(NOT InactiveCsc(Ed1)),
  count = 0;
  loop(Trays$(eff.l(Trays) > epsi AND TrayCasc(Ed1,Trays)),
    count = count + 1;
    PrevAct(Ed1,Trays,Ordr)$(ord(Ordr) eq count) = yes;
  );

);

display InitCascSize, CascSize, PrevAct;

* Reset Top Tray connectivity
loop(Ed1$(NOT InactiveCsc(Ed1)),
  loop(Trays$(TrayCasc(Ed1,Trays) AND TopTrays(Ed1,Trays)),
    OutTrV(Trays,Str)$(TopTrays(Ed1,Trays)) = no;
    OutTrV(Trays,TrayVapStr)$(ord(Trays) eq ord(TrayVapStr)) = yes;

    InTrL(Trays,Str)$(TopTrays(Ed1,Trays)) = no;
    InTrL(Trays,TrayLiqStr)$(ord(Trays) eq ord(TrayLiqStr) - 1) = yes;

* Copy stream properties/variable values from old top tray to internal stream
    loop(Str$OutTrV(Trays,Str),
      loop(Str2$OutVEd1(Ed1,Str2),
        F.l(Str) = F.l(Str2);
        Fc.l(Str,Comp) = Fc.l(Str2,Comp);
        Xc.l(Str,Comp) = Xc.l(Str2,Comp);
        Tscaled.l(Str) = Tscaled.l(Str2);
        P.l(Str) = P.l(Str2);

        ZEOS.l(Str) = ZEOS.l(Str2);
        aEOS.l(Str,j) = aEOS.l(Str2,j);
        bmEOS.l(Str) = bmEOS.l(Str2);
        amEOS.l(Str) = amEOS.l(Str2);
        bbEOS.l(Str) = bbEOS.l(Str2);
        aaEOS.l(Str) = aaEOS.l(Str2);
      );
    );

* Copy stream properties/variable values from old top tray to internal stream
    loop(Str$InTrL(Trays,Str),
      loop(Str2$InLEd1(Ed1,Str2),
        F.l(Str) = F.l(Str2);
        Fc.l(Str,Comp) = Fc.l(Str2,Comp);
        Xc.l(Str,Comp) = Xc.l(Str2,Comp);
        Tscaled.l(Str) = Tscaled.l(Str2);
        P.l(Str) = P.l(Str2);

        ZEOS.l(Str) = ZEOS.l(Str2);
        aEOS.l(Str,j) = aEOS.l(Str2,j);
        bmEOS.l(Str) = bmEOS.l(Str2);
        amEOS.l(Str) = amEOS.l(Str2);
        bbEOS.l(Str) = bbEOS.l(Str2);
        aaEOS.l(Str) = aaEOS.l(Str2);
      );
    );

  );
);

count = 0;

* Empty ActvT
ActvT(Trays) = no;

* Empty InitActive
InitActive(Trays) = no;

* Populate TrayCasc, BotTrays and TopTrays
loop(Ed1$(NOT InactiveCsc(Ed1)),
  ActvT(Trays)$(ord(Trays) > count AND ord(Trays) < count + CascSize(Ed1) + 1) = yes;
  TopTrays(Ed1,Trays) = no;
  TopTrays(Ed1,Trays)$(ord(Trays) eq count + CascSize(Ed1)) = yes;
  InitActive(Trays)$(ord(Trays) ge count + 1 AND ord(Trays) le count + InitCascSize(Ed1)) = yes;
  count = MaxStagesPerCsc + count;
);

* Establish inlet/outlet stream to top tray
loop(Ed1$(NOT InactiveCsc(Ed1)),
* Specify streams for top trays
  OutTrV(Trays,Str)$(TopTrays(Ed1,Trays)) = no;
  OutTrV(Trays,Str)$(TopTrays(Ed1,Trays) AND OutVEd1(Ed1,Str)) = yes;

  InTrL(Trays,Str)$(TopTrays(Ed1,Trays)) = no;
  InTrL(Trays,Str)$(TopTrays(Ed1,Trays) AND InLEd1(Ed1,Str)) = yes;
);

***** Temporarily Remove TrayStr from thermo calcs *****

FlashLiq(Str)$(TrayStr(Str) AND (TrayLiqStr(Str) OR TrayLiqEquil(Str))) = no;
FlashVap(Str)$(TrayStr(Str) AND (TrayVapStr(Str) OR TrayVapEquil(Str))) = no;

fEOSStr(TrayStr) = no;
RealStr(TrayStr) = no;

HCalc(Str)$TrayStr(Str) = no;
PhiCalc(Str)$TrayStr(Str) = no;

***** Repopulate TrayStr *****

* Empty TrayStr
TrayStr(Str) = no;

* Populate TrayStr
loop(ActvT,
  TrayStr(Str)$(TrayLiqEquil(Str) AND ETrL(ActvT,Str)) = yes;
  TrayStr(Str)$(TrayVapEquil(Str) AND ETrV(ActvT,Str)) = yes;

  TrayStr(Str)$(TrayLiqStr(Str) AND OutTrL(ActvT,Str)) =yes;
  TrayStr(Str)$(TrayLiqStr(Str) AND InTrL(ActvT,Str)) = yes;

  TrayStr(Str)$(TrayVapStr(Str) AND OutTrV(ActvT,Str)) = yes;
  TrayStr(Str)$(TrayVapStr(Str) AND InTrV(ActvT,Str)) = yes;
);

***** Consider New Streams for Thermo Calcs *****

FlashLiq(Str)$(TrayStr(Str) AND (TrayLiqStr(Str) OR TrayLiqEquil(Str))) = yes;
FlashVap(Str)$(TrayStr(Str) AND (TrayVapStr(Str) OR TrayVapEquil(Str))) = yes;

fEOSStr(TrayStr) = yes;
RealStr(TrayStr) = yes;

HCalc(Str)$TrayStr(Str) = yes;
PhiCalc(Str)$TrayStr(Str) = yes;


***** Initialization *****

loop(Ed1$(NOT InactiveCsc(Ed1)),
* Extract temperatures, pressures and flowrates
  loop(Str$OutVEd1(Ed1,Str),
    TTop = Tscaled.l(Str);
    PTop = P.l(Str);
    FTopV(Comp) = Fc.l(Str,Comp);
  );

  loop(Str$InLEd1(Ed1,Str),
    FTopL(Comp) = Fc.l(Str,Comp);
  );

* Grab equilibrium stream data for old top trays
  loop(Trays2$TopTraysOld(Ed1,Trays2),
    loop(Trays$TopTrays(Ed1,Trays),
* Vapor equilibrium streams
$batinclude ../ProcsFiles/CopyProps.gms ETrV(Trays2,Str2) ETrV(Trays,Str)

* Liquid equilibrium streams
$batinclude ../ProcsFiles/CopyProps.gms ETrL(Trays2,Str2) ETrL(Trays,Str)
    );
  );

  count = 0;
  loop(Trays$(TrayCasc(Ed1,Trays) AND InitActive(Trays)),
* Initialize "initially active" trays (no bypass) by shifting PrevAct tray results
    count = count + 1;
    loop(Ordr$(ord(Ordr) eq count),
      loop(Trays2$PrevAct(Ed1, Trays2, Ordr),
$ontext
        TraySrc(Trays3) = no;
        TraySrc(Trays2) = yes;
        TrayDst(Trays3) = no;
        TrayDst(Trays) = yes;
        display count, TraySrc, TrayDst;
$offtext

* OutTrV, InTrV, OutTrL, InTrL, ETrV, ETrL
* Order matters. In vapor before out vapor & out liquid before in liquid.
        InitActVarComp(Fc,InTrV);
        InitActVarComp(Fc,OutTrV);
        InitActVarComp(Fc,OutTrL);
        InitActVarComp(Fc,InTrL);
        InitActVarComp(Fc,ETrV);
        InitActVarComp(Fc,ETrL);

        InitActVarComp(Xc,InTrV);
        InitActVarComp(Xc,OutTrV);
        InitActVarComp(Xc,OutTrL);
        InitActVarComp(Xc,InTrL);
        InitActVarComp(Xc,ETrV);
        InitActVarComp(Xc,ETrL);

        InitActVar(Tscaled,InTrV);
        InitActVar(Tscaled,OutTrV);
        InitActVar(Tscaled,OutTrL);
        InitActVar(Tscaled,InTrL);
        InitActVar(Tscaled,ETrV);
        InitActVar(Tscaled,ETrL);

        InitActVar(P,InTrV);
        InitActVar(P,OutTrV);
        InitActVar(P,OutTrL);
        InitActVar(P,InTrL);
        InitActVar(P,ETrV);
        InitActVar(P,ETrL);

        InitActVar(F,InTrV);
        InitActVar(F,OutTrV);
        InitActVar(F,OutTrL);
        InitActVar(F,InTrL);
        InitActVar(F,ETrV);
        InitActVar(F,ETrL);

        InitActVar(ZEOS,OutTrV);
        InitActVar(ZEOS,InTrV);
        InitActVar(ZEOS,OutTrL);
        InitActVar(ZEOS,InTrL);
        InitActVar(ZEOS,ETrV);
        InitActVar(ZEOS,ETrL);

        eff.l(Trays) = 1;
      );
    );
  );

  loop(Trays$(TrayCasc(Ed1,Trays) AND NOT InitActive(Trays)),
* Initialize "initially inactive" trays (complete bypass - 0% efficiency)
$ontext
    Fc.l(Str,Comp)$OutTrV(Trays,Str) = FTopV(Comp);
    Fc.l(Str,Comp)$OutTrL(Trays,Str) = FTopL(Comp);

    Tscaled.l(Str)$OutTrV(Trays,Str) = TTop;
    Tscaled.l(Str)$OutTrL(Trays,Str) = TTop;

    P.l(Str)$OutTrV(Trays,Str) = PTop;
    P.l(Str)$OutTrL(Trays,Str) = PTop;
$offtext

* Vapor outlet streams
* $batinclude ../ProcsFiles/CopyProps.gms OutVEd1(Ed1,Str2) "(OutTrV(Trays,Str) OR ETrV(Trays,Str))"
$batinclude ../ProcsFiles/CopyProps.gms OutVEd1(Ed1,Str2) OutTrV(Trays,Str)

* Liquid outlet streams
$batinclude ../ProcsFiles/CopyProps.gms InLEd1(Ed1,Str2) OutTrL(Trays,Str)

    loop(Trays2$TopTrays(Ed1,Trays2),
* Problem with this implementation... ETrV, ETrL may already be overwritten

* Vapor equilibrium streams
$batinclude ../ProcsFiles/CopyProps.gms ETrV(Trays2,Str2) ETrV(Trays,Str)

* Liquid equilibrium streams
$batinclude ../ProcsFiles/CopyProps.gms ETrL(Trays2,Str2) ETrL(Trays,Str)

    );

    sL.l(Str)$(OutTrL(Trays,Str) OR ETrL(Trays,Str)) = 0;
    sV.l(Str)$(OutTrV(Trays,Str) OR ETrV(Trays,Str)) = 0;

    eff.l(Trays) = 0;
    if(TopTrays(Ed1,Trays),
      eff.fx(Trays) = 0;
    );

*    F.l(Str)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = SUM(Comp, Fc.l(Str,Comp));
*    Xc.l(Str,Comp)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = Fc.l(Str,Comp) / max(F.l(Str),1E-4);

$ontext
    Xc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Xc.l(Str2,Comp));
    Xc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Xc.l(Str2,Comp));

    Tscaled.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Tscaled.l(Str2));
    P.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), P.l(Str2));

    Tscaled.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Tscaled.l(Str2));
    P.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), P.l(Str2));

    Fc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Fc.l(Str2,Comp));
    Fc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Fc.l(Str2,Comp));

    F.l(Str)$(ETrV(Trays,Str) OR ETrL(Trays,Str)) = SUM(Comp, Fc.l(Str,Comp));
$offtext
  );
);

$batinclude ../InitFiles/Init_CEOS_Basic4.gms "(fEOSStr(Str) AND ((RealStr(Str) AND NOT InactiveStr(Str)) OR ActShdStr(Str)))"
$include ../InitFiles/Init_CEOS_Adv.gms

Kt.l(Trays,Comp)$ActvT(Trays) = SUM(Str2$ETrL(Trays,Str2), phiEOS.l(Str2,Comp))
         /max(SUM(Str$ETrV(Trays,Str), phiEOS.l(Str,Comp)),1E-5);

display InitActive, eff.l;
display TrayCasc, ActvT, AllocT, OutTrL, OutTrV, InTrL, InTrV, BotTrays, TopTrays, TrayStr;
display ETrL, ETrV;
