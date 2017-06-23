*Tscaled.l(Str) = T.l(Str)/Tref;

aEOS.l(Str,j) = Power(1 + fw(j)*(1-sqrt(Tscaled.l(Str)*Tref/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
bmEOS.l(Str) = max( SUM(j, Xc.l(Str,j)*bEOS(j)), epsi) ;
amEOS.l(Str) = max( SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2))*(1-kEOS(j,j2)))), epsi);

*amEOS.l(Str) = SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)
*    *sqrt(aEOS.l(Str,j)*aEOS.l(Str,j2)+epsi*(Power(aEOS.l(Str,j),2) + Power(aEOS.l(Str,j2),2)))*(1-kEOS(j,j2)))) ;

*bbEOS.l(Str) = max( bmEOS.l(Str)*P.l(Str)/(8.314E-2*Tscaled.l(Str)*Tref), epsi) ;
*aaEOS.l(Str) = amEOS.l(Str)*P.l(Str)/Power(8.314E-2*Tscaled.l(Str)*Tref,2) ;

bbEOS.l(Str) = max( bmEOS.l(Str)*P.l(Str)/max(8.314E-2*Tscaled.l(Str)*Tref,epsi), epsi) ;
aaEOS.l(Str) = amEOS.l(Str)*P.l(Str)/max(Power(8.314E-2*Tscaled.l(Str)*Tref,2),epsi*epsi) ;


HIG.l(Str) = 1E-3*SUM(J, Xc.l(Str,j)*(Tref**5 *CpIG('5',J)/5*(Tscaled.l(Str)**5 - 1) + Tref**4 *CpIG('4',J)/4*(Tscaled.l(Str)**4 - 1)
         + Tref**3 *CpIG('3',J)/3*(Tscaled.l(Str)**3 - 1) + Tref**2 *CpIG('2',J)/2*(Tscaled.l(Str)**2 - 1) + Tref*CpIG('1',J)*(Tscaled.l(Str) - 1)));

SIG.l(fEOSStr) = SUM(J,  Xc.l(fEOSStr, J)*(Tref**3 *CpIG('4',J)/3*(Tscaled.l(fEOSStr)**3 - 1)
         + Tref**2 *CpIG('3',J)/2*(Tscaled.l(fEOSStr)**2 - 1) + Tref*CpIG('2',J)*(Tscaled.l(fEOSStr) - 1)
         + CpIG('1',J)*log(Tscaled.l(fEOSStr)))) ;

***** Intermediate Variables *****
dadT.l(Str)$fEOSStr(Str) = -(R/2)*sqrt(omegaA/(Tscaled.l(Str)*Tref))
         *SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*(1-kEOS(j,j2))*
         (fw(j2)*sqrt(aEOS.l(Str,j)*Tc(j2)/Pc(j2)) + fw(j)*sqrt(aEOS.l(Str,j2)*Tc(j)/Pc(j)))));

delta.l(Str,j)$fEOSStr(Str) = 2*sqrt(aEOS.l(Str,j))*Sum(j2, Xc.l(Str,j2)*sqrt(aEOS.l(Str,j2))*(1 - kEOS(j, j2)))/amEOS.l(Str);

bRatio.l(Str,j)$fEOSStr(Str) = (Tc(j)/Pc(j))/max(Sum(j2, Xc.l(Str, j2)*Tc(j2)/Pc(j2)),epsi);

Inter1.l(Str)$fEOSStr(Str) = max(ZEOS.l(Str) - bbEOS.l(Str), Inter1.lo(Str));

Inter2.l(Str)$fEOSStr(Str) = max((ZEOS.l(Str) - bbEOS.l(Str))/max(ZEOS.l(Str), 1E-8), 1E-8);

Inter3.l(Str)$fEOSStr(Str) = max((2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS + Inter0))/
         (2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS - Inter0)), Inter3.lo(Str));

***** Departure Functions *****

H.l(Str)$(fEOSStr(Str) AND HCalc(Str)) = HIG.l(Str) + 0.1*( Tscaled.l(Str)*Tref*dadT.l(Str) -  amEOS.l(Str))*log(Inter3.l(Str))/(bmEOS.l(Str)*Inter0)
         + 1E-3*Rsi*Tscaled.l(Str)* Tref*(ZEOS.l(Str)-1);

S.l(Str)$(fEOSStr(Str) AND SCalc(Str)) = SIG.l(Str)
         + 100*(R*log(Inter2.l(Str)) + R*log(ZEOS.l(Str)*Pref/P.l(Str))
         + 1/(bmEOS.l(Str)*Inter0)*dadT.l(Str)*log(Inter3.l(Str)));

*phiExpPre.l(Str,j)$(fEOSStr(Str) AND PhiCalc(Str)) = bRatio.l(Str,j)*(ZEOS.l(Str)-1)-log(Inter1.l(Str))
*         + aaEOS.l(Str)/(bbEOS.l(Str)*Inter0)*(bRatio.l(Str,j) - delta.l(Str,j))*log(Inter3.l(Str));

*display Inter1.lo, Inter3.lo, Inter1.l, Inter3.l;

phiEOS.l(Str,j)$(fEOSStr(Str) AND PhiCalc(Str)) = exp(bRatio.l(Str,j)*(ZEOS.l(Str)-1)-log(Inter1.l(Str))
         + aaEOS.l(Str)/(bbEOS.l(Str)*Inter0)*(bRatio.l(Str,j) - delta.l(Str,j))*log(Inter3.l(Str))) ;

* New part added
K.l(ThrmE,j) = min(SUM(Str2$OutLGnrlE(ThrmE, Str2), phiEOS.l(Str2,j))
         /max(SUM(Str$OutVGnrlE(ThrmE,Str), phiEOS.l(Str,j)), 1E-10), 1E4);
