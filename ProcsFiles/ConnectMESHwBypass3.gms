*$offOrder
***** Establish Connectivity/Populate Sets *****

Scalar
  MaxStagesPerCsc        Maximium number of stages per cascade             ;


MaxStagesPerCsc = floor( card(Trays) / max( card(Ed1) - card(InactiveCsc), 1.0) );

Parameter
  InitCascSize(Ed1)      Initial cascade size from Edmister model solution / /;

Set
  Ordr                    Used in AdjustMESH... for shifting               /1*100/
  AllocT(Trays)           All trays allocated to a cascade                 / /
  PrevAct(Ed1,Trays,Ordr) Trays that were active in the previous solution  / / ;
*  TraySrc(Trays)
*  TrayDst(Trays)
*  StrSrc(Str)
*  StrDst(Str2);

Alias(Trays,Trays2), (Trays,Trays3);

InitCascSize(Ed1) = round(StgEd1.l(Ed1));

*loop(Ed1,
*  if(StgEd1.m(Ed1) lt 0,
*    CascSize(Ed1) = InitCascSize(Ed1) + 10;
*  else
*    CascSize(Ed1) = InitCascSize(Ed1) + 20;
*  );
*);

CascSize(Ed1) = InitCascSize(Ed1) + PlusStages;

*** Declare important macros ***
$macro InitActVar(var,setname) var.l(Str)$setname(Trays,Str) = SUM(Str2$setname(Trays2,Str2), var.l(Str2) )
$macro InitActVarComp(var,setname) var.l(Str,Comp)$setname(Trays,Str) = SUM(Str2$setname(Trays2,Str2), var.l(Str2,Comp) )

*** Systematically populate sets ***

count = 0;

* Populate TrayCasc, BotTrays and TopTrays
loop(Ed1$(NOT InactiveCsc(Ed1)),
  ActvT(Trays)$(ord(Trays) > count AND ord(Trays) < count + CascSize(Ed1) + 1) = yes;
  TrayCasc(Ed1,Trays)$(ord(Trays) > count AND ord(Trays) < count + MaxStagesPerCsc + 1) = yes;
  BotTrays(Ed1,Trays)$(ord(Trays) eq count + 1) = yes;
  TopTrays(Ed1,Trays)$(ord(Trays) eq count + CascSize(Ed1)) = yes;
  InitActive(Trays)$(ord(Trays) ge count + 1 AND ord(Trays) le count + InitCascSize(Ed1)) = yes;
  count = MaxStagesPerCsc + count;
);

* Populate AllocT
AllocT(Trays)$(ord(Trays) < count + 1) = yes;

* Establish inlet/outlet stream to tray pairs
loop(Ed1$(NOT InactiveCsc(Ed1)),
  loop(Trays$(TrayCasc(Ed1,Trays) AND NOT BotTrays(Ed1,Trays)),
    OutTrL(Trays,TrayLiqStr)$(ord(Trays) eq ord(TrayLiqStr)) = yes;
    InTrV(Trays,TrayVapStr)$(ord(Trays) eq ord(TrayVapStr)+1) = yes;
  );

  loop(Trays$(TrayCasc(Ed1,Trays) AND NOT TopTrays(Ed1,Trays)),
    OutTrV(Trays,TrayVapStr)$(ord(Trays) eq ord(TrayVapStr)) = yes;
    InTrL(Trays,TrayLiqStr)$(ord(Trays) eq ord(TrayLiqStr) - 1) = yes;
  );

* Specify streams for bottom trays
  OutTrL(Trays,Str)$(BotTrays(Ed1,Trays) AND OutLEd1(Ed1,Str)) = yes;
  InTrV(Trays,Str)$(BotTrays(Ed1,Trays) AND InVEd1(Ed1,Str)) = yes;

* Specify streams for top trays
  OutTrV(Trays,Str)$(TopTrays(Ed1,Trays) AND OutVEd1(Ed1,Str)) = yes;
  InTrL(Trays,Str)$(TopTrays(Ed1,Trays) AND InLEd1(Ed1,Str)) = yes;

  loop(Trays$TrayCasc(Ed1,Trays),
* Bypass, etc streams don't have special conditions for the top or bottom trays

    ETrL(Trays,TrayLiqEquil)$(ord(Trays) eq ord(TrayLiqEquil)) = yes;
    ETrV(Trays,TrayVapEquil)$(ord(Trays) eq ord(TrayVapEquil)) = yes;
  );
);

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

display TrayCasc, ActvT, AllocT, OutTrL, OutTrV, InTrL, InTrV, BotTrays, TopTrays, TrayStr;
display ETrL, ETrV;

* Remove exit casacde stream from InThrmEPres(Flsh,Str)
* This prevents a redundent pressure equation
*loop(Ed1,
*  loop(Flsh,
*    InThrmEPres(Flsh,Str)$(OutLEd1(Ed1,Str)) = no;
*  );
*);

