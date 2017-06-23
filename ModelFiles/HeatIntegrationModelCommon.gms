* These are common files for the updated and revised heat integration models.

*Alias (PStr, PStr2);

Positive Variables
  Tin(GnrlE)                      Inlet temperature
  Tout(GnrlE)                     Outlet temperature

  TPcnd(Str,HIZone)              Pinch candidate temperatures
  QsZ(HIZone)                    Steam (hot) utility load for each zone
  QwZ(HIZone)                    Water (cold) utility load for each zone
  Qs                             Total hot utility
  Qw                             Total cold utility;

Equations
  EqQwZ(HIZone)                  Calculate cooling water demand from energy balance + Qs
  EqQsTotal                      Calculate total hot utility
  EqQwTotal                      Calculate total cold utility
  EqMaxDT_cool(HtEx)             Specify max change in temperature for cooling HXs
  EqMaxDT_heat(HtEx)             Specify max change in temperature for heating HXs
  EqMaxDT_cond(TCond)            Specify max change in temperature for condensers
  EqMaxDT_reb(PReb)              Specify max change in temperature for reboilers
  ;

EqQwZ(HIZone)..
  QwZ(HIZone) =e= QsZ(HIZone) + SUM(GnrlE$(NOT InactiveGnrlE(GnrlE) AND HIMap(GnrlE, HIZone)), Qout(GnrlE) - Qin(GnrlE));

EqQsTotal..
  Qs =e= SUM(HIZone, QsZ(HIZone));

EqQwTotal..
  Qw =e= SUM(HIZone, QwZ(HIZone));

EqMaxDT_cool(HtEx)$CoolHtEx(HtEx)..
  SUM(Str$InOneGnrlE(HtEx,Str), Tscaled(Str)) - SUM(Str2$OutVHtEx(HtEx,Str2), Tscaled(Str2)) =l= MaxDT/Tref;

EqMaxDT_heat(HtEx)$HeatHtEx(HtEx)..
  SUM(Str2$OutVHtEx(HtEx,Str2), Tscaled(Str2)) - SUM(Str$InOneGnrlE(HtEx,Str), Tscaled(Str)) =l= MaxDT/Tref;

EqMaxDT_cond(TCond)..
  SUM(Str$InOneGnrlE(TCond,Str), Tscaled(Str)) - SUM(Str2$OutTCond(TCond,Str2), Tscaled(Str2)) =l= MaxDT/Tref;

EqMaxDT_reb(PReb)..
  SUM(Str2$OutVPReb(PReb,Str2), Tscaled(Str2)) - SUM(Str$InOneGnrlE(PReb,Str), Tscaled(Str)) =l= MaxDT/Tref;
