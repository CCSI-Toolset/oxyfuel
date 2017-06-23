********************************
***** Init_CEOS_Basic4.gms *****
********************************

* This file initializes ZEOS using analytical solutions. Note: Use $batinclude
* Input arguments
*
* Input 1: conditionals (index Str) to select streams to initialize.
* Typical input: (fEOSStr(Str) AND ((RealStr(Str) AND NOT InactiveStr(Str)) OR ActShdStr(Str)))
*
* Input 2: To be added?


$ontext
loop(Str,
  if(SUM(Comp, Xc.l(Str,Comp)) < 0.99 OR SUM(Comp, Xc.l(Str,Comp)) > 1.01,
      Xc.lo(Str,Comp) = 0;
      Xc.up(Str,Comp) = 1;
      Xc.l(Str,Comp) = FeedFlow(Comp) / SUM(j, FeedFlow(j));
  );
);
$offtext

aEOS.l(Str,j)$%1 = Power(1 + fw(j)*(1-sqrt(Tscaled.l(Str)*Tref/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
bmEOS.l(Str)$%1 = max( SUM(j, Xc.l(Str,j)*bEOS(j)), epsi) ;
amEOS.l(Str)$%1 = max( SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2))*(1-kEOS(j,j2)))), epsi);

*amEOS.l(Str) = SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)
*    *sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2)+epsi*(Power(aEOS.l(Str,j),2) + Power(aEOS.l(Str,j2),2)))*(1-kEOS(j,j2)))) ;

bbEOS.l(Str)$%1 = max( bmEOS.l(Str)*P.l(Str)/(8.314E-2*Tscaled.l(Str)*Tref), epsi) ;
aaEOS.l(Str)$%1 = amEOS.l(Str)*P.l(Str)/Power(8.314E-2*Tscaled.l(Str)*Tref,2) ;


*ZEOS.l(Str) = 0.001;
*ZEOS.l(LiqStr) = bbEOS.l(LiqStr)+0.0001;
*ZEOS.l(VapStr) = 0.8;

display FlashVap, FlashLiq, F.l;

loop(Str$%1,
  StrNow(Str2) = no;
  StrNow(Str) = yes;

  if(DebugMode > 2,
    display StrNow;
  );
  aZ = -(1 + bbEOS.l(Str) - uEOS*bbEOS.l(Str));
  bZ = aaEOS.l(Str) + wEOS*Power(bbEOS.l(Str),2) - uEOS*bbEOS.l(Str) - uEOS*Power(bbEOS.l(Str),2);
  cZ = - aaEOS.l(Str)*bbEOS.l(Str) - wEOS*Power(bbEOS.l(Str),2) - wEOS*Power(bbEOS.l(Str),3);
  Qzi = (aZ*aZ - 3*bZ)/9;
  Rzi = (2*Power(aZ,3) - 9*aZ*bZ + 27*cZ)/54;
  Mzi = Power(Rzi,2) - Power(Qzi,3);

  if(DebugMode > 2,
    display Mzi;
  );

  if(Mzi <= 0,
    ThetaZi = arccos(Rzi/rPower(Qzi,3/2));
    Zi('1') = -(2*sqrt(Qzi)*cos(ThetaZi/3) ) - aZ/3;
    Zi('2') = -(2*sqrt(Qzi)*cos((ThetaZi + 2*pi)/3) ) - aZ/3;
    Zi('3') = -(2*sqrt(Qzi)*cos((ThetaZi - 2*pi)/3) ) - aZ/3;

    if(DebugMode > 2,
      display Zi;
    );

    ZEOS.l(Str)$(FlashVap(Str) OR VapStr(Str)) = smax(RootsZi$(3*Power(Zi(RootsZi),2) - 2*(1 + (1- uEOS)*bbEOS.l(Str))*Zi(RootsZi)
         + ( aaEOS.l(Str) - uEOS*bbEOS.l(Str) + (wEOS - uEOS)*Power(bbEOS.l(Str),2) ) > 0  ), Zi(RootsZi));

* AND Zi(RootsZi) > bbEOS.l(Str)

    ZEOS.l(Str)$(FlashLiq(Str) OR LiqStr(Str)) = smin(RootsZi$(3*Power(Zi(RootsZi),2) - 2*(1 + (1- uEOS)*bbEOS.l(Str))*Zi(RootsZi)
         + ( aaEOS.l(Str) - uEOS*bbEOS.l(Str) + (wEOS - uEOS)*Power(bbEOS.l(Str),2) ) > 0 ), Zi(RootsZi));

*    ZEOS.l(Str)$(FlashLiq(Str) OR LiqStr(Str)) = smin(RootsZi$(Zi(RootsZi) > 0), Zi(RootsZi));

* Check if the calculations are close to the three vs one root cutoff
* If so, use F as a guide
    if(abs(Mzi) < 1E-8 AND F.l(Str) > 1E-4,

* Classified as a vapor, but identified a liquid root
      if( (FlashVap(Str) OR VapStr(Str)) AND 6*ZEOS.l(Str) - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) < epsi,
        ZEOS.l(Str) = 0.8;

* Classified as a liquid, but identified a vapor root
      elseif (FlashLiq(Str) OR LiqStr(Str)) AND 6*ZEOS.l(Str) - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) > -epsi,
        ZEOS.l(Str) = bbEOS.l(Str) + 0.001;
      );
    );

  else
    Szi = -sign(Rzi)*rPower(abs(Rzi) + sqrt(Mzi), 1/3);
    if(Szi eq 0,
      Tzi = 0;
    else
      Tzi = Qzi/Szi;
    );
    Zi('1') = Szi + Tzi - aZ/3;
    Zi('2') = 0;
    Zi('3') = 0;
    ZEOS.l(Str) = Zi('1');
*    if( Zi('1') < 0.4 AND ( ( F.l(Str) > 0 AND (FlashVap(Str) OR VapStr(Str)) ) OR ( RelaxThermo(Str) AND FlashLiq(Str) ) ),
*      ZEOS.l(Str) = 0.8;
*    elseif Zi('1') > 0.6 AND ((F.l(Str) > 0 AND (FlashLiq(Str) OR LiqStr(Str))) OR (RelaxThermo(Str) AND FlashVap(Str) )),
*      ZEOS.l(Str) = bbEOS.l(Str) + 0.001;
*    );

    if(F.l(Str) > 1E-3 AND 6*Zi('1') - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) < -epsi AND (FlashVap(Str) OR VapStr(Str)),
      ZEOS.l(Str) = 0.8;
    elseif 6*Zi('1') - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) < epsi AND RelaxThermo(Str) AND FlashLiq(Str),
      ZEOS.l(Str) = 0.8;
    elseif F.l(Str) > 1E-3 AND 6*Zi('1') - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) > epsi AND (FlashLiq(Str) OR LiqStr(Str)),
      ZEOS.l(Str) = bbEOS.l(Str) + 0.001;
    elseif 6*Zi('1') - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) > -epsi AND RelaxThermo(Str) AND FlashVap(Str),
      ZEOS.l(Str) = bbEOS.l(Str) + 0.001;
    );

    if(DebugMode > 2,
      display Zi;
    );
* display Szi, Tzi, Qzi, Rzi, Mzi;

* Adjust sL and sV if F = 0 and phase should be relaxed
    if(F.l(Str) < 1E-4 AND FlashVap(Str) AND 6*ZEOS.l(Str) - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) < 0,
      sV.l(Str) = 10;
    elseif F.l(Str) < 1E-4 AND FlashLiq(Str) AND 6*ZEOS.l(Str) - 2*(1 + (1 - uEOS)*bbEOS.l(Str)) > 0,
      sL.l(Str) = 10;
    );
  );
);

if(DebugMode > 1,
  display T.l, P.l, sL.l, sV.l, ZEOS.l;
);

V.l(Str)$%1 = ZEOS.l(Str)*R*Tscaled.l(Str)*Tref/max(P.l(Str),1E-6);
