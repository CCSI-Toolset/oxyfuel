* This file initializes ZEOS using analytical solutions. Note: Use $batinclude
* Input arguments
*
* Input 1: conditionals (index Str) to select streams to initialize.


HIG.l(Str)$%1 = 1E-3*SUM(J, Xc.l(Str,j)*(Tref**5 *CpIG('5',J)/5*(Tscaled.l(Str)**5 - 1) + Tref**4 *CpIG('4',J)/4*(Tscaled.l(Str)**4 - 1)
         + Tref**3 *CpIG('3',J)/3*(Tscaled.l(Str)**3 - 1) + Tref**2 *CpIG('2',J)/2*(Tscaled.l(Str)**2 - 1) + Tref*CpIG('1',J)*(Tscaled.l(Str) - 1)));

SIG.l(Str)$(%1 AND fEOSStr(Str)) = SUM(J,  Xc.l(Str, J)*(Tref**3 *CpIG('4',J)/3*(Tscaled.l(Str)**3 - 1)
         + Tref**2 *CpIG('3',J)/2*(Tscaled.l(Str)**2 - 1) + Tref*CpIG('2',J)*(Tscaled.l(Str) - 1)
         + CpIG('1',J)*log(Tscaled.l(Str)))) ;

***** Intermediate Variables *****
dadT.l(Str)$(%1 AND fEOSStr(Str)) = -(R/2)*sqrt(omegaA/(Tscaled.l(Str)*Tref))
         *SUM(j, SUM(j2, Xc.l(Str,j)*Xc.l(Str,j2)*(1-kEOS(j,j2))*
         (fw(j2)*sqrt(aEOS.l(Str,j)*Tc(j2)/Pc(j2)) + fw(j)*sqrt(aEOS.l(Str,j2)*Tc(j)/Pc(j)))));

delta.l(Str,j)$(%1 AND fEOSStr(Str)) = 2*sqrt(aEOS.l(Str,j))*Sum(j2, Xc.l(Str,j2)*sqrt(aEOS.l(Str,j2))*(1 - kEOS(j, j2)))/amEOS.l(Str);

bRatio.l(Str,j)$(%1 AND fEOSStr(Str)) = (Tc(j)/Pc(j))/max(Sum(j2, Xc.l(Str, j2)*Tc(j2)/Pc(j2)),epsi);

Inter1.l(Str)$(%1 AND fEOSStr(Str)) = max(ZEOS.l(Str) - bbEOS.l(Str), Inter1.lo(Str));

Inter2.l(Str)$(%1 AND fEOSStr(Str)) = max((ZEOS.l(Str) - bbEOS.l(Str))/max(ZEOS.l(Str), 1E-8), 1E-8);

Inter3.l(Str)$(%1 AND fEOSStr(Str)) = max((2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS + Inter0))/
         (2*ZEOS.l(Str) + bbEOS.l(Str)*(uEOS - Inter0)), Inter3.lo(Str));

***** Departure Functions *****

H.l(Str)$(%1 AND fEOSStr(Str) AND HCalc(Str)) = HIG.l(Str) + 0.1*( Tscaled.l(Str)*Tref*dadT.l(Str) -  amEOS.l(Str))*log(Inter3.l(Str))/(bmEOS.l(Str)*Inter0)
         + 1E-3*Rsi*Tscaled.l(Str)* Tref*(ZEOS.l(Str)-1);

S.l(Str)$(%1 AND fEOSStr(Str) AND SCalc(Str)) = SIG.l(Str)
         + 100*(R*log(Inter2.l(Str)) + R*log(ZEOS.l(Str)*Pref/P.l(Str))
         + 1/(bmEOS.l(Str)*Inter0)*dadT.l(Str)*log(Inter3.l(Str)));

phiEOS.l(Str,j)$(%1 AND fEOSStr(Str) AND PhiCalc(Str)) = exp(bRatio.l(Str,j)*(ZEOS.l(Str)-1)-log(Inter1.l(Str))
         + aaEOS.l(Str)/(bbEOS.l(Str)*Inter0)*(bRatio.l(Str,j) - delta.l(Str,j))*log(Inter3.l(Str))) ;

loop(ThrmE$(NOT InactiveGnrlE(ThrmE)),
  K.l(ThrmE,j) = min(SUM(Str2$OutLGnrlE(ThrmE, Str2), phiEOS.l(Str2,j))
         /max(SUM(Str$OutVGnrlE(ThrmE,Str), phiEOS.l(Str,j)), 1E-10), 1E4);
);
