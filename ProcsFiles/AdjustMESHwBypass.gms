NumTrays.l(Ed1) = SUM(Trays$TrayCasc(Ed1,Trays), eff.l(Trays));

InitCascSize(Ed1) = round(NumTrays.l(Ed1));
CascSize(Ed1) = min(InitCascSize(Ed1) + PlusStages, MaxStagesPerCsc);

display InitCascSize, CascSize;

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
    TTop = T.l(Str);
    PTop = P.l(Str);
    FTopV(Comp) = Fc.l(Str,Comp);
  );

  loop(Str$InLEd1(Ed1,Str),
    FTopL(Comp) = Fc.l(Str,Comp);
  );

  loop(Trays$(TrayCasc(Ed1,Trays) AND NOT InitActive(Trays)),
* Initialize "initially inactive" trays (complete bypass - 0% efficiency)
    Fc.l(Str,Comp)$OutTrV(Trays,Str) = FTopV(Comp);
    Fc.l(Str,Comp)$OutTrL(Trays,Str) = FTopL(Comp);

    T.l(Str)$OutTrV(Trays,Str) = TTop;
    T.l(Str)$OutTrL(Trays,Str) = TTop;

    P.l(Str)$OutTrV(Trays,Str) = PTop;
    P.l(Str)$OutTrL(Trays,Str) = PTop;

    loop(Str2$OutVEd1(Ed1,Str2),
      loop(Str$(OutTrV(Trays,Str) OR ETrV(Trays,Str)),
        ZEOS.l(Str) = ZEOS.l(Str2);
        aEOS.l(Str,j) = aEOS.l(Str2,j);
        bmEOS.l(Str) = bmEOS.l(Str2);
        amEOS.l(Str) = amEOS.l(Str2);
        bbEOS.l(Str) = bbEOS.l(Str2);
        aaEOS.l(Str) = aaEOS.l(Str2);
      );
    );

    loop(Str2$OutLEd1(Ed1,Str2),
      loop(Str$(OutTrL(Trays,Str) OR ETrL(Trays,Str)),
        ZEOS.l(Str) = ZEOS.l(Str2);
        aEOS.l(Str,j) = aEOS.l(Str2,j);
        bmEOS.l(Str) = bmEOS.l(Str2);
        amEOS.l(Str) = amEOS.l(Str2);
        bbEOS.l(Str) = bbEOS.l(Str2);
        aaEOS.l(Str) = aaEOS.l(Str2);
      );
    );

    eff.l(Trays) = 0;

    F.l(Str)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = SUM(Comp, Fc.l(Str,Comp));
    Xc.l(Str,Comp)$(OutTrV(Trays,Str) OR OutTrL(Trays,Str)) = Fc.l(Str,Comp) / max(F.l(Str),1E-4);

    Xc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Xc.l(Str2,Comp));
    Xc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Xc.l(Str2,Comp));

    T.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), T.l(Str2));
    P.l(Str)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), P.l(Str2));

    T.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), T.l(Str2));
    P.l(Str)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), P.l(Str2));

    Fc.l(Str,Comp)$ETrV(Trays,Str) = SUM(Str2$OutTrV(Trays, Str2), Fc.l(Str2,Comp));
    Fc.l(Str,Comp)$ETrL(Trays,Str) = SUM(Str2$OutTrL(Trays, Str2), Fc.l(Str2,Comp));

    F.l(Str)$(ETrV(Trays,Str) OR ETrL(Trays,Str)) = SUM(Comp, Fc.l(Str,Comp));
  );
);

$include ../InitFiles/Init_CEOS_Adv.gms

Kt.l(Trays,Comp)$ActvT(Trays) = SUM(Str2$ETrL(Trays,Str2), phiEOS.l(Str2,Comp))
         /max(SUM(Str$ETrV(Trays,Str), phiEOS.l(Str,Comp)),1E-5);

display InitActive;
display TrayCasc, ActvT, AllocT, OutTrL, OutTrV, InTrL, InTrV, BotTrays, TopTrays, TrayStr;
display ETrL, ETrV;
