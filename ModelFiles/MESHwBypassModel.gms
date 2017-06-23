Sets
  TrayVapBypass(Str)     Streams used for vapor bypass           /Vb1*Vb250/
  TrayLiqBypass(Str)     Streams used for liquid bypass          /Lb1*Lb250/
  TrayVapIn(Str)         Streams into equil section (vapor)      /Vin1*Vin250/
  TrayLiqIn(Str)         Streams into equil section (liquid)     /Lin1*Lin250/
  TrayVapEquil(Str)      Streams out of equil section (vapor)    /Ve1*Ve250/
  TrayLiqEquil(Str)      Streams out of equil section (liquid)   /Le1*Le250/

  ByTrV(Trays,Str)       Bypass stream-tray pairing (vapor)      / /
  ByTrL(Trays,Str)       Bypass stream-tray pairing (liquid)     / /
  InSecV(Trays,Str)      Section inlet stream-tray pairing (vap) / /
  InSecL(Trays,Str)      Section inlet stream-tray pairing (liq) / /
  OutSecV(Trays,Str)     Section outlet stream-tray pairing (vap)/ /
  OutSecL(Trays,Str)     Section outlet stream-tray pairing (liq)/ /

  ByVSptr(Sptr)          Bypass splitters                       /Ssv1*Ssv250/
  ByLSptr(Sptr)          Bypass splitters                       /Ssl1*Ssl250/
  ByVMxr(Mxr)            Bypass mixers                          /Mv1*Mv250/
  ByLMxr(Mxr)            Bypass mixers                          /Ml1*Ml250/

  InitActive(Trays)      Trays that are initially active;

***** Equations *****

Alias(Trays, Trays2);

Positive Variables
  effV(Trays)
  effL(Trays)
  eff(Trays);

Equations
  EqCalcTrayEffV(Trays,Str,Str2)
  EqCalcTrayEffL(Trays,Str,Str2)
  EqCalcOverallEff(Trays)
  EqOrderEfficiencies(Ed1, Trays, Trays2);

EqCalcTrayEffV(Trays,Str,Str2)$(InSecV(Trays,Str) AND InTrV(Trays,Str2))..
  F(Str) =e= effV(Trays)*F(Str2);

EqCalcTrayEffL(Trays,Str,Str2)$(InSecL(Trays,Str) AND InTrL(Trays,Str2))..
  F(Str) =e= effL(Trays)*F(Str2);

EqCalcOverallEff(Trays)..
  eff(Trays) =e= (effL(Trays) + effV(Trays))/2;

EqOrderEfficiencies(Ed1, Trays, Trays2)$(TrayCasc(Ed1,Trays) AND TrayCasc(Ed1,Trays2) AND ord(Trays) gt 1 AND ord(Trays2) eq ord(Trays) + 1)..
  eff(Trays) =g= eff(Trays2);

Model BypassEffEqns /EqCalcTrayEffV, EqCalcTrayEffL, EqCalcOverallEff, EqOrderEfficiencies/;