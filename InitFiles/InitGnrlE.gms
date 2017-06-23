Q(GnrlE) = SUM(Str$(OutLGnrlE(GnrlE, Str) OR OutVGnrlE(GnrlE, Str)), F.l(Str)*H.l(Str))
  - SUM(Str$InGnrlE(GnrlE, Str), F.l(Str)*H.l(Str));

loop(GnrlE,
  if(Q(GnrlE) > 0,
    Qin.l(GnrlE) = Q(GnrlE);
    Qout.l(GnrlE) = 0;
  else
    Qin.l(GnrlE) = 0;
    Qout.l(GnrlE) = -Q(GnrlE);
  );
);
