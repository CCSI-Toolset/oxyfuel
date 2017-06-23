Parameter chckP(Str);

chckP(Str) = Pb.l(Str) - Pd.l(Str);

loop(Str,
  if(SUM(Comp, Xc.l(Str,Comp)) < 0.99 OR SUM(Comp, Xc.l(Str,Comp)) > 1.01,
      Xc.lo(Str,Comp) = 0;
      Xc.up(Str,Comp) = 1;
      Xc.l(Str,Comp) = FeedFlow(Comp) / SUM(j, FeedFlow(j));
  );
);


aEOS.l(Str,j) = Power(1 + fw(j)*(1-sqrt(Tscaled.l(Str)*Tref/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
bmEOS.l(Str) = max( SUM(j, Xc.l(Str,j)*bEOS(j)), epsi) ;
amEOS.l(Str) = max( SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2))*(1-kEOS(j,j2)))), epsi);

*amEOS.l(Str) = SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)
*    *sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2)+epsi*(Power(aEOS.l(Str,j),2) + Power(aEOS.l(Str,j2),2)))*(1-kEOS(j,j2)))) ;

bbEOS.l(Str) = max( bmEOS.l(Str)*P.l(Str)/(8.314E-2*Tscaled.l(Str)*Tref), epsi) ;
aaEOS.l(Str) = amEOS.l(Str)*P.l(Str)/Power(8.314E-2*Tscaled.l(Str)*Tref,2) ;

ZEOS.l(Str) = 0.001;
ZEOS.l(LiqStr) = bbEOS.l(LiqStr)+0.0001;
ZEOS.l(VapStr) = 0.8;

display FlashVap, FlashLiq;

loop(Str$FlashVap(Str),
*  if(sV.l(Str) <= 0.01,
  if(P.l(Str) <= Pd.l(Str) + epsi OR chckP(Str) < 0,
    ZEOS.l(Str) = 0.8;
  else
    ZEOS.l(Str) = bbEOS.l(Str)+0.0001;
  );
);

loop(Str$FlashLiq(Str),
*  if(sL.l(Str) <= 0.01,
  if(P.l(Str) >= Pb.l(Str) - epsi OR chckP(Str) < 0,
    ZEOS.l(Str) = bbEOS.l(Str)+0.0001;
  else
    ZEOS.l(Str) = 0.8;
  );
);

display T.l, P.l, Pb.l, Pd.l, sL.l, sV.l, ZEOS.l;

V.l(Str) = ZEOS.l(Str)*R*T.l(Str)/max(P.l(Str),1E-6);
