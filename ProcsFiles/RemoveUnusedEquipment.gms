*************************************
***** RemoveUnusedEquipment.gms *****
*************************************

* This file is used to prune unused equipment, helping reduce unneeded equations
* and degeneracy. It is intended to be run after optimization is completed but before
* the solution is refined with additional optimizaiton runs.

* For every GnrlE if there is zero feed flow, prune the unit and its outlet
* streams.
loop(GnrlE,
  if(SUM(Str$InGnrlE(GnrlE,Str), F.l(Str)) eq 0,
    InactiveGnrlE(GnrlE) = yes;
    InactiveStr(Str)$(OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = yes;
    Qin.fx(GnrlE) = 0;
    Qout.fx(GnrlE) = 0;
*    F.fx(Str) = 0;
  );
);

* For every pruned streams fix the flow at zero and remove the streams to RealStr.
F.fx(Str)$InactiveStr(Str) = 0;
Fc.fx(Str,Comp)$InactiveStr(Str) = 0;
RealStr(Str)$InactiveStr(Str) = no;

* Don't consider inactive streams for thermo calculations
* fEOSStr(Str)$InactiveStr(Str) = no;

* Remove pruned streams from heat integration consideration.
*PStr(Str)$InactiveStr(Str) = no;
*InPStr(Str)$InactiveStr(Str) = no;
*OutPStr(Str)$InactiveStr(Str) = no;

display InactiveStr, InactiveGnrlE, PStr, InPStr, OutPStr;
