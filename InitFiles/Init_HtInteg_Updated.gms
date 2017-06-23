* Heating Units (aka cold streams)
Tin.l(GnrlE)$HEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), Tref*Tscaled.l(Str));
*Tin.l(GnrlE)$HEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), T.l(Str)) + alpha;
*Tin.l(GnrlE)$HEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), sum(Str2$OutLGnrlE(GnrlE,Str2), T.l(Str2) - smmax(T.l(Str2) - T.l(Str) + alpha) + alpha ));

Tout.l(GnrlE)$HEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), sum(Str2$OutLGnrlE(GnrlE,Str2), smmax(Tref*Tscaled.l(Str2) - Tref*Tscaled.l(Str) + alpha) + Tref*Tscaled.l(Str) - alpha ));
*Tout.l(GnrlE)$HEqp(GnrlE) = sum(Str2$OutLGnrlE(GnrlE,Str2), T.l(Str2));

* Cooling Units (aka hot streams)
Tin.l(GnrlE)$CEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), sum(Str2$OutLGnrlE(GnrlE,Str2), smmax(Tref*Tscaled.l(Str) - Tref*Tscaled.l(Str2) + alpha) + Tref*Tscaled.l(Str2) - alpha ));
*Tin.l(GnrlE)$CEqp(GnrlE) =  sum(Str$InOneGnrlE(GnrlE,Str), T.l(Str));
*Tin.l(GnrlE)$CEqp(GnrlE) = sum(Str$InOneGnrlE(GnrlE,Str), sum(Str2$OutLGnrlE(GnrlE,Str2), T.l(Str2) + smmax(T.l(Str) - T.l(Str2) + alpha) - alpha ));

Tout.l(GnrlE)$CEqp(GnrlE) = sum(Str2$OutLGnrlE(GnrlE,Str2), Tref*Tscaled.l(Str2));

* Initialize FCp
FCp.l(GnrlE)$(HEqp(GnrlE)) = Qin.l(GnrlE)/(max(Tout.l(GnrlE) - Tin.l(GnrlE), 1E-8));
FCp.l(GnrlE)$(CEqp(GnrlE)) = Qout.l(GnrlE)/(max(Tin.l(GnrlE) - Tout.l(GnrlE), 1E-8));


* Initialize the pinch points
loop(HIZone,
  loop(GnrlE$(CEqp(GnrlE) AND HIMap(GnrlE, HIZone)),
    TPcnd.l(Str,HIZone)$(InOneGnrlE(GnrlE,Str) AND PStr(Str, HIZone)) = Tin.l(GnrlE);
    TPcnd.l(Str,HIZone)$(OutLGnrlE(GnrlE,Str) AND PStr(Str, HIZone)) = Tout.l(GnrlE);
  );
);

loop(HIZone,
  loop(GnrlE$(HEqp(GnrlE) AND HIMap(GnrlE, HIZone)),
    TPcnd.l(Str,HIZone)$(InOneGnrlE(GnrlE,Str) AND PStr(Str,HIZone)) = Tin.l(GnrlE) + HRAT(HIZone);
    TPcnd.l(Str,HIZone)$(OutLGnrlE(GnrlE,Str) AND PStr(Str,HIZone)) = Tout.l(GnrlE) + HRAT(HIZone);
  );
);

loop(HIZone,
loop(Str$PStr(Str,HIZone),
*  QAhU.l(Str) = SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE)), FCp.l(GnrlE)*smmax(Tin.l(GnrlE) - TPcnd.l(Str)) - FCp.l(GnrlE)*smmax( Tout.l(GnrlE) - TPcnd.l(Str) ) );
*  QAcU.l(Str) = SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE)), FCp.l(GnrlE)*smmax( Tout.l(GnrlE) - TPcnd.l(Str) + HRAT) - FCp.l(GnrlE)*smmax(Tin.l(GnrlE) - TPcnd.l(Str) + HRAT ) );

QAhU.l(Str,HIZone)$(PStr(Str, HIZONE) AND NOT InactiveStr(Str)) = SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,Str) AND NOT OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)),
                         FCp.l(GnrlE)*smmax(Tin.l(GnrlE) - TPcnd.l(Str, HIZone)) - FCp.l(GnrlE)*smmax( Tout.l(GnrlE) - TPcnd.l(Str, HIZone) ) )
                   + SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)), Qout.l(GnrlE));

QAcU.l(Str,HIZone)$(PStr(Str, HIZone) AND NOT InactiveStr(Str)) = SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,Str) AND NOT OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)),
                         FCp.l(GnrlE)*smmax( Tout.l(GnrlE) - TPcnd.l(Str, HIZone) + HRAT(HIZone)) - FCp.l(GnrlE)*smmax(Tin.l(GnrlE) - TPcnd.l(Str, HIZone) + HRAT(HIZone) ) )
                   + SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND InOneGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)), Qin.l(GnrlE));

);
);

*display Qin.l, Qout.l, QAhU.l, QAcU.l, FCp.l;

QsZ.l(HIZone) = max(smax(Str$PStr(Str,HIZone), QAcU.l(Str,HIZone) - QAhU.l(Str, HIZone)), 0);

QwZ.l(HIZone) = max(QsZ.l(HIZone) + SUM(GnrlE$(NOT InactiveGnrlE(GnrlE) AND HIMap(GnrlE, HIZone)), Qout.l(GnrlE) - Qin.l(GnrlE)),0);

Qs.l = SUM(HIZone, QsZ.l(HIZone));
Qw.l = SUM(HIZone, QwZ.l(HIZone));
