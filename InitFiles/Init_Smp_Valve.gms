loop(Str,
  PinV.l(Valve)$InValve(Valve,Str) = P.l(Str);
  TinV.l(Valve)$InValve(Valve,Str) = Tscaled.l(Str)*Tref;
);

loop(Str,
  PoutV.l(Valve)$OutLValve(Valve,Str) = P.l(Str);
  ToutV.l(Valve)$OutLValve(Valve,Str) = Tscaled.l(Str)*Tref;
);

dPv.l(Valve) = PoutV.l(Valve) - PinV.l(Valve);

Pv.l(Valve) = 0.5*(PinV.l(Valve) + PoutV.l(Valve));

dTv.l(Valve) = ToutV.l(Valve) - TinV.l(Valve);

Tv.l(Valve) = 0.5*(TinV.l(Valve) + ToutV.l(Valve));

XcV.l(Valve,Comp) = SUM(Str$InValve(Valve,Str), Fc.l(Str,Comp))/
                         max(SUM(Comp2, SUM(Str$InValve(Valve,Str), Fc.l(Str,Comp2))), 1E-10);

CpRef.l(Valve) = SUM(j, XcV.l(Valve,j)*(CpIG('5',j)*Tv.l(Valve)**4 + CpIG('4',j)*Tv.l(Valve)**3 + CpIG('3',j)*Tv.l(Valve)**2 +
                         CpIG('2',j)*Tv.l(Valve) + CpIG('1',j)));

alphaV.l(Valve) = (2.1417 - 0.8413*Pv.l(Valve) - 0.0356*(Tv.l(Valve) - 273.15) - 0.0094*(Tv.l(Valve) - 273.15)*Pv.l(Valve))/1000;

dTv.l(Valve) = Rsi*( Tv.l(Valve) )*( alphaV.l(Valve)*Tv.l(Valve) - 1 )*dPv.l(Valve)
         /max( CpRef.l(Valve)*Pv.l(Valve), 1E-10 );
