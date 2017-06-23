**********************************
***** CommonThermoModels.gms *****
**********************************

* This file contains equations required for both the Simple and CEOS
* thermodynamic models.

*----------------------------------------------------*
* Ideal Gas Thermodynamics
*----------------------------------------------------*

* equation for molar enthalpy in units of kJ/mol I think
Equation
EqHIG(Str)       Ideal gas enthalpy;

EqHIG(Str)$(fEOSStr(Str) AND HCalc(Str) AND RealStr(Str))..
  HIG(Str) =e= 1E-3*SUM(J, Xc(Str,j)*(Tref**5 *CpIG('5',J)/5*(Tscaled(Str)**5 - 1) + Tref**4 *CpIG('4',J)/4*(Tscaled(Str)**4 - 1)
 + Tref**3 *CpIG('3',J)/3*(Tscaled(Str)**3 - 1) + Tref**2 *CpIG('2',J)/2*(Tscaled(Str)**2 - 1) + Tref*CpIG('1',J)*(Tscaled(Str) - 1)));

* equation for molar entropy in units of J/(mol K) I think
Equation
EqSIG(Str)       Ideal gas entropy;

EqSIG(Str)$(fEOSStr(Str) AND SCalc(Str)  AND RealStr(Str))..
  SIG(Str) =e= SUM(J, Xc(Str,j)*(
                 Tref**3*CpIG('4',J)/3*(Tscaled(Str)**3 - 1)
                 + Tref**2 *CpIG('3',J)/2*(Tscaled(Str)**2 - 1)
                 + Tref*CpIG('2',J)*(Tscaled(Str) - 1)
                 + CpIG('1',J)*log(Tscaled(Str))
                 )) ;

Equations
  EqSlackL(Str,ThrmE)    Relaxation for equilibrium for vanishing liquid streams
  EqSlackV(Str,ThrmE)    Relaxation for equilibrium for vanishing vapor streams
  EqSlack(Str,Str2,ThrmE) Equality constraint version;

EqSlackL(Str,ThrmE)$(OutLGnrlE(ThrmE,Str) AND NOT InactiveGnrlE(ThrmE))..
  -sL(Str) =l= beta(ThrmE) - 1;

EqSlackV(Str,ThrmE)$(OutVGnrlE(ThrmE,Str) AND NOT InactiveGnrlE(ThrmE))..
   beta(ThrmE) - 1 =l= sV(Str);

EqSlack(Str,Str2,ThrmE)$(OutLGnrlE(ThrmE,Str) AND OutVGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  beta(ThrmE) =e= 1 - sL(Str) + sV(Str);
