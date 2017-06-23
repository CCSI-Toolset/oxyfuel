*******************************
***** CEOSThermoModel.gms *****
*******************************

* This file contains the cubic EOS thermo model equations. Details regarding the
* root selection formulation are discuss in
*
* "An equation-oriented approach for handling thermodynamics based
* on cubic equation of state in process optimization" by Kamath, Biegler
* and Grossmann (2010)
*
* Departure function formulas are taken from
* "The Properties of Gases and Liquids" by Reid, Prausnitz and Sherwood, 4th Ed

Parameters
  EpsilonZ       Smallest value for dfdZ in CEOS                 /1E-6/
  EpsilonA       Smallest value for dfdA in CEOS                 /1E-6/

Variables
 bmEOS(Str)        bm for mixture for CEOS
 phiEOS(Str,j)     Fugacity coefficient
* phiExpPre(Str,j)  pre-exponential term for fugacity calculations
 ZEOS(Str)         Compressibility factor for CEOS
 bbEOS(Str)        big B for CEOS
 aEOS(Str,j)       Final a for CEOS
 amEOS(Str)        am for mixture for CEOS
 aaEOS(Str)        big A for CEOS
 dadT(Str)         T times partial a over partial T
 delta(Str,j)      Intermeidate term for fugacity departure function
 bRatio(Str,j)     b(j) over b ratio for fugacity departure function
 Inter1(Str)       Intermediate term for departure function (Z - B)
 Inter2(Str)       Intermediate term for departure function (Z + B)
 Inter2(Str)       Intermediate term for departure function (Z - B) over Z
 Inter4(Str)       Intermediate term for departure function (Z + B) over Z
 Inter3(Str)       Intermediate term for departure function (2Z + B(u + Inter0)) over (2Z + B(u - Inter0))
;

*----------------------------------------------------*
* Basic Equations for CEOS
*----------------------------------------------------*
Equations
EqaEOSf(Str,j)   Calculate final a for CEOS
EqbmEOSf(Str)    Calculate bm for mixture
EqamEOSf(Str)    Calculate am for mixture
EqbbEOSf(Str)    Calculate big B
EqaaEOSf(Str)    Calculate big A
EqZEOSf(Str)     Cubic equation of state
EqZf(Str)        Compressibility factor definition
EqEOSf(Str)      First derivative constraint
EqVEOSf(Str)     Second derivative contraint for vapor streams (non vanishing)
EqLEOSf(Str)     Second derivative contraint for liqud streams (non vanishing)
EqVEOSflash(Str) Second derivative contraint for vapor streams (vanishing)
EqLEOSflash(Str) Second derivative contraint for liquid streams (vanishing)
EqLEOSflash2(Str)        Helps prevent divide by zero - obsolete
EqlinZLEOSf(Str)         Helps prevent divide by zero - obsolete;

$macro cubicCond(junk) fEOSStr(Str) AND ((RealStr(Str) AND NOT InactiveStr(Str)) OR ActShdStr(Str))

EqaEOSf(Str,j)$(cubicCond(0))..
  aEOS(Str,j) =E= Power(1 + fw(j)*(1-sqrt(Tscaled(Str)*Tref/Tc(j))),2)*omegaA*Power(8.314E-2*Tc(j),2)/Pc(j) ;
* Unit of a : bar (m3/kmol)^2
EqbmEOSf(Str)$(cubicCond(0))..
  bmEOS(Str) =E= SUM(j, Xc(Str,j)*bEOS(j)) ;
* Unit of bEOS: m^3/kmol

EqamEOSf(Str)$(cubicCond(0))..
* used Prof. Biegler's modification for case when a goes to zero
*  amEOS(Str) =E= SUM(j, SUM(j2, Xc(Str,j)*Xc(Str,j2)
*    *sqrt(aEOS(Str,j)*aEOS(Str,j2)+epsi*(Power(aEOS(Str,j),2) + Power(aEOS(Str,j2),2)))*(1-kEOS(j,j2)))) ;
amEOS(Str) =E= SUM(j, SUM(j2, Xc(Str,j)*Xc(Str,j2)*sqrt(aEOS(Str,j)*aEOS(Str,j2))*(1-kEOS(j,j2))));

EqbbEOSf(Str)$(cubicCond(0))..
  bbEOS(Str) =E= bmEOS(Str)*P(Str)/(8.314E-2*Tscaled(Str)*Tref) ;

EqaaEOSf(Str)$(cubicCond(0))..
  aaEOS(Str) =E= amEOS(Str)*P(Str)/Power(8.314E-2*Tscaled(Str)*Tref,2) ;

EqZEOSf(Str)$(cubicCond(0))..
  (Power(ZEOS(Str),3) - (1 + bbEOS(Str) - uEOS*bbEOS(Str))*Power(ZEOS(Str),2)
         + (aaEOS(Str) + (wEOS - uEOS)*Power(bbEOS(Str),2) - uEOS*bbEOS(Str))*ZEOS(Str)
         - aaEOS(Str)*bbEOS(Str) - wEOS*Power(bbEOS(Str),2) - wEOS*Power(bbEOS(Str),3))
         =E= 0 ;

EqZf(Str)$(cubicCond(0) AND VCalc(Str))..
  R*Tscaled(Str)*Tref*ZEOS(Str) =E= P(Str)*V(Str) ;

EqEOSf(Str)$(cubicCond(0))..
  3*Power(ZEOS(Str),2) - 2*(1 + (1- uEOS)*bbEOS(Str))*ZEOS(Str)
         + ( aaEOS(Str) - uEOS*bbEOS(Str) + (wEOS - uEOS)*Power(bbEOS(Str),2) ) =G= EpsilonZ ;

EqVEOSf(Str)$(VapStr(Str) AND cubicCond(0))..
  6*ZEOS(Str) - 2*(1 + (1 - uEOS)*bbEOS(Str))  =G= 0 ;

EqLEOSf(Str)$(LiqStr(Str) AND cubicCond(0))..
  6*ZEOS(Str) - 2*(1 + (1 - uEOS)*bbEOS(Str))  =L= 0 ;

EqVEOSflash(Str)$(FlashVap(Str) AND cubicCond(0))..
  6*ZEOS(Str) - 2*(1 + (1 - uEOS)*bbEOS(Str))  =G= -bigM*sV(Str) ;

EqLEOSflash(Str)$(FlashLiq(Str) AND cubicCond(0))..
  6*ZEOS(Str) - 2*(1 + (1 - uEOS)*bbEOS(Str))  =L= bigM*sL(Str) ;

EqLEOSflash2(Str)$(FlashLiq(Str) AND cubicCond(0))..
  ZEOS(Str) - bbEOS(Str) =G= 1.0E-5 ;

EqlinZLEOSf(Str)$(LiqStr(Str) AND cubicCond(0))..
  ZEOS(Str) - bbEOS(Str) =G= 1.0E-5 ;

Model BasicCEOSEqns /EqaEOSf, EqbmEOSf, EqamEOSf, EqbbEOSf, EqaaEOSf, EqZEOSf,
         EqZf, EqEOSf, EqLEOSf, EqVEOSf, EqVEOSflash, EqLEOSflash,
         EqSlackL, EqSlackV/;
*         EqSlack/;

sL.fx(Str)$(NOT FlashLiq(Str)) = 0;
sV.fx(Str)$(NOT FlashVap(Str)) = 0;

*----------------------------------------------------*
* Phase Stability
*----------------------------------------------------*

Equations
  EqPhaseStabLiq(Str,Str2,Str3)
  EqPhaseStabFlashLiq(Str,Str2,Str3)
  EqPhaseStabVap(Str,Str2,Str3)
  EqPhaseStabFlashVap(Str,Str2,Str3);

EqPhaseStabLiq(Str,VapShd,LiqShd)$(ComMap(Str,VapShd,LiqShd) AND PhaseStability(Str)
         AND LiqStr(Str))..
  Tscaled(Str) =l= Tscaled(LiqShd);

EqPhaseStabFlashLiq(Str,VapShd,LiqShd)$(ComMap(Str,VapShd,LiqShd) AND PhaseStability(Str)
         AND FlashLiq(Str))..
  Tscaled(Str) =l= Tscaled(LiqShd) + sL(Str)/Tref;

EqPhaseStabVap(Str,VapShd,LiqShd)$(ComMap(Str,VapShd,LiqShd) AND PhaseStability(Str)
         AND VapStr(Str))..
  Tscaled(Str) =g= Tscaled(VapShd);

EqPhaseStabFlashVap(Str,VapShd,LiqShd)$(ComMap(Str,VapShd,LiqShd) AND PhaseStability(Str)
         AND FlashVap(Str))..
  Tscaled(Str) =g= Tscaled(VapShd) - sV(Str)/Tref;

*----------------------------------------------------*
* Departure Function Equations for CEOS
*----------------------------------------------------*

* Intermediate equations used for departure function calculations
Equations
  EqdadT(Str)
  Eqdelta(Str,j)
  EqbRatio(Str,j)
  EqInter1(Str)
  EqInter2(Str)
  EqInter3(Str);

EqdadT(Str)$(cubicCond(0))..
  dadT(Str)*sqrt(Tscaled(Str)*Tref) =e= -(R/2)*sqrt(omegaA)*SUM(j, SUM(j2, Xc(Str,j)*Xc(Str,j2)*(1-kEOS(j,j2))*
         (fw(j2)*sqrt(aEOS(Str,j)*Tc(j2)/Pc(j2)) + fw(j)*sqrt(aEOS(Str,j2)*Tc(j)/Pc(j)))));
* unit of dadT : (m^3/kmol)^2 * bar /K

Eqdelta(Str,j)$(cubicCond(0))..
  delta(Str,j)*amEOS(Str) =e= 2*sqrt(aEOS(Str,j))*Sum(j2, Xc(Str,j2)*sqrt(aEOS(Str,j2))*(1 - kEOS(j, j2)));

EqbRatio(Str,j)$(cubicCond(0))..
  bRatio(Str,j)*Sum(j2, Xc(Str, j2)*Tc(j2)/Pc(j2)) =e= Tc(j)/Pc(j);

EqInter1(Str)$(cubicCond(0))..
  Inter1(Str) =e= ZEOS(Str) - bbEOS(Str);

EqInter2(Str)$(cubicCond(0))..
  Inter2(Str)*ZEOS(Str) =e= ZEOS(Str) - bbEOS(Str);

EqInter3(Str)$(cubicCond(0))..
  Inter3(Str)*(2*ZEOS(Str) + bbEOS(Str)*(uEOS - Inter0))
         =e= (2*ZEOS(Str) + bbEOS(Str)*(uEOS + Inter0));

Inter1.lo(Str) = 10**(EpsilonA);
Inter2.lo(Str) = epsi;
Inter3.lo(Str) = 0.5;

Model IntermEOSEqns /EqdadT, Eqdelta, EqbRatio, EqInter1, EqInter2, EqInter3/;

* equation for molar enthalpy in units of kJ/mol I think
Equations
EqHEOS(Str)      Ethalpy Departure function;

EqHEOS(Str)$(cubicCond(0) AND HCalc(Str))..
H(Str)*bmEOS(Str) =e= HIG(Str)*bmEOS(Str) + 0.1*( Tscaled(Str)*Tref*dadT(Str)-amEOS(Str) )*log(Inter3(Str))/Inter0
                      + 1E-3*Rsi*Tscaled(Str)*Tref*(ZEOS(Str)-1)*bmEOS(Str);

* Units for H are kJ/mol
* Units of amEOS/bmEOS are m^3*bar/kmol... multiplied by 100 yields J/mol
* m^3*bar = 100,000 J = 100 kJ
* m^3*bar/kmol = 100,000 J/kmol = 100 J /mol = 0.1 kJ/mol

Model HEOSEqns /EqHIG, EqHEOS /;

* equation for molar entropy in units of J/(mol K) I think
Equations
EqSEOS(Str)      Entropy departure function;

EqSEOS(Str)$(cubicCond(0) AND SCalc(Str))..
  S(Str) =e= SIG(Str) + 100*(
         R*log(Inter2(Str)) + R*log(ZEOS(Str)*Pref/P(Str))
         + 1/(bmEOS(Str)*Inter0)*dadT(Str)*log(Inter3(Str)) );

Model SEOSEqns /EqSEOS, EqSIG/;

* equations for fugacity
Equations
EqphiEOSf(Str,j)         Fugacity equation
EqphiExp(Str,j)          Calculate Fugacity intermediate (phiExpPre)
EqPhiK(ThrmE,Str,Str2,j) Calculate equilibrium constant (K) from fugacities;

EqphiEOSf(Str,j)$(cubicCond(0) AND PhiCalc(Str))..
  (phiEOS(Str,j)*Inter1(Str))**(bbEOS(Str)*Inter0)
         =E= Inter3(Str)**(aaEOS(Str)*(bRatio(Str,j) - delta(Str,j)) )*
                 exp(bRatio(Str,j)*( ZEOS(Str)-1 )*bbEOS(Str)*Inter0);

*  phiEOS(Str,j) =E= exp(bRatio(Str,j)*(ZEOS(Str)-1)-log(Inter1(Str))
*         + aaEOS(Str)/(bbEOS(Str)*Inter0)*(bRatio(Str,j) - delta(Str,j))
*  *log(Inter3(Str))) ;


*EqphiExp(Str,j)$(fEOSStr(Str) AND PhiCalc(Str))..
*  phiExpPre(Str,j) =e= bRatio(Str,j)*(ZEOS(Str)-1)-log(Inter1(Str))
*         + aaEOS(Str)/(bbEOS(Str)*Inter0)*(bRatio(Str,j) - delta(Str,j))*log(Inter3(Str));

***** Bound fugacities *****
*phiExpPre.up(Str,j) = log(1E5);
*phiExpPre.lo(Str,j) = log(1E-5);

EqPhiK(ThrmE,Str,Str2,j)$(OutVGnrlE(ThrmE,Str) AND OutLGnrlE(ThrmE, Str2) AND NOT InactiveGnrlE(ThrmE))..
  K(ThrmE,j)*phiEOS(Str,j) =e= phiEOS(Str2,j);

phiEOS.lo(Str,j) = epsi;

Model CubicEOSEqns /BasicCEOSEqns, IntermEOSEqns, HEOSEqns, SEOSEqns, EqphiEOSf, EqPhiK/;
* EqphiExp
