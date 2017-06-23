*******************************
***** EntropyBalModel.gms *****
*******************************

* This file contains "lite" entropy balances - in other words only entropy
* balances for the valves and cascades. Typically this file is not used.

* Originally it was noticed that including an entropy balance for valves helped
* convergence. With refinements of the initialization procedure, this constraint
* has been removed.
*
* The cascade entropy balance should be unneccissary, but is availabe in case
* the Edmister equations are not trusted. A better approach is to calculate the
* entropy generation in the cascades post-optimization. (In entropy generation
* is negative there is a problem.)

*** "Lite" Entropy Balances ***

Equations
*  EqEntpBalGnrlE(GnrlE)
  EqEntpBalCsc(Ed1)
  EqEntpBalValve(Valve)
;

*EqEntpBalGnrlE(GnrlE)$(Cmpr(GnrlE))..
*Pump(GnrlE) or
*  SUM(Str$InGnrlE(GnrlE,Str), F(Str)*S(Str)) =l= SUM(Str2$(OutVGnrlE(GnrlE,Str2) or OutLGnrlE(GnrlE,Str2)), F(Str2)*S(Str2));

EqEntpBalCsc(Ed1)..
  SUM(Str$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)), F.l(Str)*S.l(Str)) =l= SUM(Str2$(OutVEd1(Ed1,Str2) OR OutLEd1(Ed1,Str2)), F.l(Str2)*S.l(Str2));

EqEntpBalValve(Valve)..
  SUM(Str$InValve(Valve,Str), F(Str)*S(Str)) =l= SUM(Str2$(OutVValve(Valve,Str2) or OutLValve(Valve,Str2)), F(Str2)*S(Str2));
*EqEntpBalCsc
*Model EntropyBalsLite /EqEntpBalGnrlE/;

*Model EntropyBals /EqEntpBalThrmE, EqEntpBalQTCond/;
