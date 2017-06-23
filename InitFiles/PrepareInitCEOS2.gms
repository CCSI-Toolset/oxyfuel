********************************
****** PrepareInitCEOS.gms *****
********************************

* This file handles the intermediate processing between solving the simplified
* thermo model ASU and initialize the CEOS variables.


***** Enable Shadow Streams *****
* The shadow streams (used for BP/DP calc) are not included in the simple thermo
* ASU model. Vapor pressure correlations are instead used to ensure BP or DP.
* The shadow streams much be initialized and added into the optimization problem.

* Initialize shadow stream compositions, temperatures and pressures
loop(Str$(DewPoint(Str) OR BubPoint(Str)),
  loop(Str2,
    loop(Str3$((ComMap(Str,Str2,Str3) AND BubPoint(Str)) OR (ComMap(Str,Str3,Str2) AND DewPoint(Str))),
      Xc.l(Str2,Comp) = Xc.l(Str,Comp);
      Fc.l(Str2,Comp) = Xc.l(Str,Comp);
      Tscaled.l(Str2) = Tscaled.l(Str);
      P.l(Str2) = P.l(Str);
    );
  );
);

* Create advanced only streams

$ontext
Sets
  AdvOnly(Str);

AdvOnly(VapStrAdvOnly) = yes;
AdvOnly(LiqStrAdvOnly) = yes;

display AdvOnly;

* Consider AdvOnly streams for thermo calculations

fEOSStr(AdvOnly) = yes;
LiqStr(LiqStrAdvOnly) = yes;
VapStr(VapStrAdvOnly) = yes;

loop(Str$PhaseStability(Str),
  loop(LiqShd,
    loop(VapShd$ComMap(Str,VapShd,LiqShd),
      AdvOnly(LiqShd) = yes;
      AdvOnly(VapShd) = yes;
    );
  );
);

fEOSStr(AdvOnly) = yes;

LiqStr(Str)$(ActShdStr(Str) AND LiqShd(Str)) = yes;
VapStr(Str)$(ActShdStr(Str) AND VapShd(Str)) = yes;
$offtext
*phiEOS.l(Str,Comp) = 1;
*
*T.fx(FeedStr) = FeedT;

*T.fx(FeedStr) = 300;
