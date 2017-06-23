*********************************
***** SimpleThermoModel.gms *****
*********************************

* This file contains the correlation based thermodynamic model used
* primarily for initialization. Phases are determined from vapor pressure,
* bubble point and dew point calculations. Stream ethanlpy is calculated from
* an ideal gas polynomial and the Watson correlation. Cooling in valves is
* approximated using a correlation for the Joule-Tompson coefficient.

Variables
   TLiqSlack(Str,Comp)           Slack used for Watson correlation - avoids negative logs
   XcV(Valve,Comp)               Composition entering a valve
   CpRef(Valve)                  Ideal gas heat capacity
   Tv(Valve)                     Average temperature in valve
   Pv(Valve)                     Average pressure in valve
   dTv(Valve)                    Temperature change in valve
   dPv(Valve)                    Pressure drop accross valve
   alphaV(Valve)                 Joule-Tompson coefficient for valve;

Equations
  EqHLiqAlt(Str)                 Calculate enthalpy for liquid streams (alternate form)
  EqHVap(Str)                    Calculate enthalpy for vapor streams
  EqHVapSurf(Str)                Surface fit using CEOS data from Aspen
  EqHLiqSurf(Str)                Surface fit using CEOS data from Aspen
  EqPBub(Str)                    Calculate bubble point for stream
  EqPDew(Str)                    Calculate dew point for stream
  EqPBubDewChck(Str)             Ensure Pd <= Pb
  EqSimpleBubPoint(Str)          Constraint specified streams at bubble point
  EqSimpleDewPoint(Str)          Constraint specified streams at dew point
  EqSimpleL(Str)                 Phase constraint for liquid streams (non-vanishing)
  EqSimpleRelaxL(Str)            Phase constraint for liquid streams (vanishing)
  EqSimpleV(Str)                 Phase constraint for vapor streams (non-vanishing)
  EqSimpleRelaxV(Str)            Phase constraint for vapor streams (vanishing)
  EqPvap(Str,Comp)               Vapor pressure calculations
  EqSimpleK(ThrmE,Str,Comp)      Calculate phase equilibrium constant
  EqThrmEDew(ThrmE,Str,Str2)     sL calculation - not used
  EqThrmEBub(ThrmE,Str,Str2)     sV calculation - not used
  EqXcV(Valve,Comp)              Calculate XcV - valve correlation model
  EqCpRef(Valve)                 Calculate CpRef - valve correlation model
  EqTv(Valve)                    Calculate Tv - valve correlation model
  EqPv(Valve)                    Calculate Pv - valve correlation model
  EqdT(Valve)                    Calculate dT - valve correlation model
  EqdP(Valve)                    Calculate dP - valve correlation model
  EqalphaV(Valve)                Calculate alphaV - valve correlation model
  EqValveTempCalc(Valve)         Joule-Tompson effect - relates dT and dP
  EqRevisedJTEffect(Valve)
  EqDoubleRevisedJTEffect(Valve)

  ;

$macro sThermCond(junk) RealStr(Str) AND NOT InactiveStr(Str)

EqHVap(Str)$((VapStr(Str) OR FlashVap(Str)) AND sThermCond(0))..
  H(Str) =e= HIG(Str);

EqHLiqAlt(Str)$((LiqStr(Str) OR FlashLiq(Str)) AND sThermCond(0))..
  H(Str) =e= HIG(Str) - SUM(Comp, Xc(Str,Comp)*1E-3*Rsi*(AntConst('B',Comp)/Power(Tscaled(Str)*Tref + AntConst('C',Comp),2) ));

EqHVapSurf(Str)$((VapStr(Str) OR FlashVap(Str)) AND sThermCond(0) AND HCalc(Str))..
  H(Str) =e= SUM(Comp, Xc(Str,Comp)*
         (HVapSurf('1',Comp) + HVapSurf('2',Comp)*P(Str) + HVapSurf('3',Comp)*Power(P(Str),2)
         + HVapSurf('4',Comp)*Tscaled(Str) + HVapSurf('5',Comp)*Power(Tscaled(Str),2)
         + HVapSurf('6',Comp)*Power(Tscaled(Str),3) + HVapSurf('7',Comp)*P(Str)*Tscaled(Str)
         + HVapSurf('8',Comp)*P(Str)*Power(Tscaled(Str),2)));

EqHLiqSurf(Str)$((LiqStr(Str) OR FlashLiq(Str)) AND sThermCond(0) AND HCalc(Str))..
  H(Str) =e= SUM(Comp, Xc(Str,Comp)*(HLiqSurf('1',Comp)
         + HLiqSurf('4',Comp)*Tscaled(Str) + HLiqSurf('5',Comp)*Power(Tscaled(Str),2)
         + HLiqSurf('6',Comp)*Power(Tscaled(Str),3) + HLiqSurf('7',Comp)*P(Str)*Tscaled(Str)
         + HLiqSurf('8',Comp)*P(Str)*Power(Tscaled(Str),2)));


EqPBub(Str)$RealStr(Str)..
    Pb(Str) =e= sum(Comp, Xc(Str,Comp)*Pvap(Str,Comp));

Pvap.lo(Str,Comp) = 0.01;

EqPDew(Str)$RealStr(Str)..
    Pd(Str)*sum(Comp, Xc(Str,Comp)/( Pvap(Str,Comp) ) ) =e= 1;
*Pd(Str)*SUM(Comp, Xc(Str,Comp)*Prod(Comp2$(ord(Comp2) <> ord(Comp)), Pvap(Str,Comp2) )) =e= Prod(Comp, Pvap(Str,Comp));

EqPBubDewChck(Str)$RealStr(Str)..
  Pd(Str) =l= Pb(Str);

EqSimpleL(Str)$(LiqStr(Str) AND sThermCond(0) AND NOT BubPoint(Str))..
P(Str) =g= Pb(Str);

EqSimpleRelaxL(Str)$(FlashLiq(Str) AND sThermCond(0))..
P(Str) =g= Pb(Str) - bigM*sL(Str);

EqSimpleV(Str)$(VapStr(Str) AND sThermCond(0) AND NOT DewPoint(Str))..
P(Str) =l= Pd(Str);

EqSimpleRelaxV(Str)$(FlashVap(Str) AND sThermCond(0))..
P(Str) =l= Pd(Str) + bigM*sV(Str);

EqSimpleBubPoint(Str)$(BubPoint(Str) AND sThermCond(0))..
  Pb(Str) =e= P(Str);

EqSimpleDewPoint(Str)$(DewPoint(Str) AND sThermCond(0))..
  Pd(Str) =e= P(Str);

EqPvap(Str, Comp)$(sThermCond(0))..
* Pvap in bar after conversion from mmHg
  Pvap(Str, Comp) =e= exp(AntConst('A',Comp) - AntConst('B',Comp)/(Tscaled(Str)*Tref + AntConst('C',Comp)))/750.061683;

EqSimpleK(ThrmE,Str,Comp)$(OutVGnrlE(ThrmE,Str) AND NOT InactiveGnrlE(ThrmE))..
  K(ThrmE,Comp)*P(Str) =e= Pvap(Str,Comp);

EqThrmEDew(ThrmE,Str,Str2)$(OutVGnrlE(ThrmE,Str) and OutLGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  P(Str) - sL(Str2) =e= Pd(Str);

EqThrmEBub(ThrmE,Str,Str2)$(OutVGnrlE(ThrmE,Str) and OutLGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  P(Str2) + sV(Str) =e= Pb(Str2);

EqXcV(Valve,Comp)$(NOT InactiveGnrlE(Valve))..
  XcV(Valve,Comp)*SUM(Comp2, SUM(Str$InValve(Valve,Str), Fc(Str,Comp2))) =e= SUM(Str$InValve(Valve,Str), Fc(Str,Comp));

EqTv(Valve)$(NOT InactiveGnrlE(Valve))..
  2*Tv(Valve) =e= TinV(Valve) + ToutV(Valve);

EqPv(Valve)$(NOT InactiveGnrlE(Valve))..
  2*Pv(Valve) =e= PinV(Valve) + PoutV(Valve);

EqCpRef(Valve)$(NOT InactiveGnrlE(Valve))..
  CpRef(Valve) =e= SUM(j, XcV(Valve,j)*(CpIG('5',j)*Tv(Valve)**4 + CpIG('4',j)*Tv(Valve)**3 + CpIG('3',j)*Tv(Valve)**2 +
                         CpIG('2',j)*Tv(Valve) + CpIG('1',j)));

EqdP(Valve)$(NOT InactiveGnrlE(Valve))..
  dPv(Valve) =e= PoutV(Valve) - PinV(Valve);

EqdT(Valve)$(NOT InactiveGnrlE(Valve))..
  dTv(Valve) =e= ToutV(Valve) - TinV(Valve);

EqalphaV(Valve)$(NOT InactiveGnrlE(Valve))..
  1000*alphaV(Valve) =e= 2.1417 - 0.8413*Pv(Valve)
         - 0.0356*(Tv(Valve) - 273.15)
         - 0.0094*(Tv(Valve) - 273.15)*Pv(Valve);

EqValveTempCalc(Valve)$(NOT InactiveGnrlE(Valve))..
  dTv(Valve)*CpRef(Valve)*Pv(Valve) =e= Rsi*( Tv(Valve) )*( alphaV(Valve)*Tv(Valve) - 1 )*dPv(Valve);
*  dTv(Valve)*0.5*(P(Str)+P(Str2)) =e= Rsi*0.042*dPv(Valve);

*EqRevisedJTEffect(Valve)$(NOT InactiveGnrlE(Valve))..
*  ToutV(Valve) - TinV(Valve) =e= JTCoeff*(PoutV(Valve) - PinV(Valve));

EqDoubleRevisedJTEffect(Valve)$(NOT InactiveGnrlE(Valve))..
  ToutV(Valve) - TinV(Valve) =e= 0.299853*exp(-0.009125*TinV(Valve))*(PoutV(Valve) - PinV(Valve));

dTv.lo(Valve) = -25;

*Model SimpleThermoValve /EqCpRef, EqXcV, EqdP, EqPv, EqdT, EqTv, EqalphaV, EqValveTempCalc/;
* EqalphaV, EqValveTempCalc
Model SimpleThermoValve /EqDoubleRevisedJTEffect/;
* EqRevisedJTEffect

*Model SimpleThermo / EqHIG, EqHLiqAlt, EqHVap, EqPBub, EqPDew, EqSimpleL, EqSimpleRelaxL,
*         EqSimpleV, EqSimpleRelaxV, EqSimpleK, EqPvap, EqSlackL, EqSlackV, SimpleThermoValve/;
* EqPBubDewChck
* EqHLiq, EqHLiqF, EqLiqSlack1, EqLiqSlack2,

Model SimpleThermo / EqHVapSurf, EqHLiqSurf, EqPBub, EqPDew, EqSimpleL, EqSimpleRelaxL,
         EqSimpleV, EqSimpleRelaxV, EqSimpleK, EqPvap, EqSlackL, EqSlackV, EqSimpleBubPoint, EqSimpleDewPoint/;
