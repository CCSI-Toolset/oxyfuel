Sets
  TrayVapEquil(Str)      Streams out of equil section (vapor)    /Ve1*Ve350/
  TrayLiqEquil(Str)      Streams out of equil section (liquid)   /Le1*Le350/

  ETrV(Trays,Str)        Equilibrium stream (vapor)      / /
  ETrL(Trays,Str)        Equilibrium stream (liquid)     / /

  InitActive(Trays)      Trays that are initially active;

***** Equations *****

Alias(Trays, Trays2);

Variables
  Kslack(Trays,Comp);

Positive Variables
  betaT(Trays);

Equations
  EqTrayKCEOS2(Trays,Comp,Str,Str2)
  EqTrayKCEOS2Alt(Trays,Comp,Str,Str2)
  EqTrayEnrgBalMxrV(Trays)
  EqTrayEnrgBalMxrL(Trays)
  EqTrayMassBalMxrV(Trays,Str,Str2,Str3,Comp)
  EqTrayMassBalMxrL(Trays,Str,Str2,Str3,Comp)
  EqTrayMassBal2(Trays,Comp)
  EqTrayEnrgBal2(Trays)
  EqTrayXcEquilib2(Trays,Comp,Str,Str2)
  EqTrayTEquilib2(Trays,Str,Str2)
  EqTrayPEquilib2(Trays,Str,Str2)
  EqOrderEfficiencies(Ed1, Trays, Trays2)
  EqNumStages(Ed1)
  EqSlackTL(Trays,Str)
  EqSlackTV(Trays,Str)
;

EqTrayKCEOS2Alt(Trays,Comp,Str,Str2)$(ETrV(Trays,Str) AND ETrL(Trays,Str2) AND ActvT(Trays))..
  Kt(Trays,Comp)*phiEOS(Str,Comp) =e= phiEOS(Str2,Comp) + Kslack(Trays,Comp);

EqTrayKCEOS2(Trays,Comp,Str,Str2)$(ETrV(Trays,Str) AND ETrL(Trays,Str2) AND ActvT(Trays))..
  Kt(Trays,Comp)*phiEOS(Str,Comp) =e= phiEOS(Str2,Comp);

EqTrayEnrgBalMxrV(Trays)$ActvT(Trays)..
  SUM(Str$OutTrV(Trays,Str), F(Str)*H(Str)) =e=
         eff(Trays)*SUM(Str2$ETrV(Trays,Str2), F(Str2)*H(Str2))
         + (1 - eff(Trays))*SUM(Str3$InTrV(Trays,Str3), F(Str3)*H(Str3));

EqTrayEnrgBalMxrL(Trays)$ActvT(Trays)..
  SUM(Str$OutTrL(Trays,Str), F(Str)*H(Str)) =e=
         eff(Trays)*SUM(Str2$ETrL(Trays,Str2), F(Str2)*H(Str2))
         + (1 - eff(Trays))*SUM(Str3$InTrL(Trays,Str3), F(Str3)*H(Str3));

EqTrayMassBalMxrV(Trays,Str,Str2,Str3,Comp)$(InTrV(Trays,Str) AND OutTrV(Trays,Str2) AND ETrV(Trays,Str3) AND ActvT(Trays))..
  Fc(Str2,Comp) =e= eff(Trays)*Fc(Str3,Comp) + (1 - eff(Trays))*Fc(Str,Comp);

EqTrayMassBalMxrL(Trays,Str,Str2,Str3,Comp)$(InTrL(Trays,Str) AND OutTrL(Trays,Str2) AND ETrL(Trays,Str3) AND ActvT(Trays))..
  Fc(Str2,Comp) =e= eff(Trays)*Fc(Str3,Comp) + (1 - eff(Trays))*Fc(Str,Comp);

EqTrayMassBal2(Trays,Comp)$ActvT(Trays)..
  SUM(Str$(InTrL(Trays,Str) OR InTrV(Trays,Str)), Fc(Str,Comp)) =e=
   SUM(Str2$(ETrL(Trays,Str2) OR ETrV(Trays,Str2)), Fc(Str2,Comp));

EqTrayEnrgBal2(Trays)$ActvT(Trays)..
  SUM(Str$(InTrL(Trays,Str) OR InTrV(Trays,Str)), F(Str)*H(Str)) =e=
    SUM(Str2$(ETrL(Trays,Str2) OR ETrV(Trays,Str2)), F(Str2)*H(Str2));

EqTrayXcEquilib2(Trays,Comp,Str,Str2)$(ETrV(Trays,Str) AND ETrL(Trays,Str2) AND ActvT(Trays))..
  Xc(Str,Comp) =e= betaT(Trays)*Kt(Trays,Comp)*Xc(Str2,Comp);

EqTrayTEquilib2(Trays,Str,Str2)$(ETrV(Trays,Str) AND ETrL(Trays,Str2) AND ActvT(Trays))..
  Tscaled(Str) =e= Tscaled(Str2);

EqTrayPEquilib2(Trays,Str,Str2)$((ETrV(Trays,Str) AND OutTrV(Trays,Str2)) OR (ETrL(Trays,Str) AND OutTrL(Trays,Str2)) AND ActvT(Trays))..
  P(Str) =e= P(Str2);

EqOrderEfficiencies(Ed1, Trays, Trays2)$(TrayCasc(Ed1,Trays) AND TrayCasc(Ed1,Trays2) AND ord(Trays) lt card(Trays) AND ord(Trays2) eq ord(Trays) + 1 AND ActvT(Trays))..
  eff(Trays) =g= eff(Trays2);

EqNumStages(Ed1)..
  NumTrays(Ed1) =e= SUM(Trays$(TrayCasc(Ed1,Trays) AND ActvT(Trays)), eff(Trays));

EqSlackTL(Trays,Str)$(ETrL(Trays,Str) AND ActvT(Trays))..
  -sL(Str) =l= betaT(Trays) - 1;

EqSlackTV(Trays,Str)$(ETrV(Trays,Str) AND ActvT(Trays))..
   betaT(Trays) - 1 =l= sV(Str);

eff.lo(Trays) = 0;
eff.up(Trays) = 1;

Model TrayEqnsWithBypass /EqTrayPEquilib, EqTrayPDrop, EqTrayConstPres,
                         EqTrayKCEOS2, EqTrayEnrgBalMxrV, EqTrayEnrgBalMxrL, EqTrayMassBalMxrV, EqTrayMassBalMxrL,
                         EqTrayMassBal2, EqTrayEnrgBal2, EqTrayXcEquilib2, EqTrayTEquilib2, EqTrayPEquilib2, EqNumStages, EqSlackTL, EqSlackTV/;

Model MassBalancesForTrayModel /EqTrayPEquilib, EqTrayPDrop, EqTrayConstPres, EqTrayMassBalMxrV, EqTrayMassBalMxrL,
                         EqTrayMassBal2, EqTrayTEquilib2, EqTrayPEquilib2/;