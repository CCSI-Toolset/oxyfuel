parameter dT /75/;

loop(HtEx$CoolHtEx(HtEx),
  T.l(Str)$(OutLHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str)) = SUM(Str2$InHtExOne(HtEx,Str2), T.l(Str2)) - dT;
  P.l(Str)$(OutLHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str)) = SUM(Str2$InHtExOne(HtEx,Str2), P.l(Str2));
  Qout.l(HtEx) = 10;
);

loop(HtEx$HeatHtEx(HtEx),
  T.l(Str)$(OutLHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str)) = SUM(Str2$InHtExOne(HtEx,Str2), T.l(Str2)) + dT;
  P.l(Str)$(OutLHtEx(HtEx,Str) OR OutVHtEx(HtEx,Str)) = SUM(Str2$InHtExOne(HtEx,Str2), P.l(Str2));
  Qin.l(HtEx) = 10;
);

Tscaled.l(Str) = T.l(Str)/Tref;
display T.l, P.l;
