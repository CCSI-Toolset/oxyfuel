**************************
***** ThrmEModel.gms *****
**************************

* This file contains general equations for the "thermodynamic
* equipment models". All ThrmE units share a similar structure:
*   - One or more feeds
*   - Only two possible exit streams: liquid and vapor
*   - Exit streams are in thermodynamic equilibrium

* To differentiate each type of ThrmE equipment, additional
* temperature and pressure relationships are required. These
* are specified in seperate model files.

*----------------------------------------------------*
*        Variables for General Thermo Equipment
*----------------------------------------------------*

Positive Variables
 Qin(GnrlE)              Heat added for each ThrmE (kJ per time)
 Qout(GnrlE)             Heat removed for each ThrmE (kJ per time);

Variables
 beta(ThrmE)             Used for thermodynamic relaxation
 K(ThrmE,Comp)           Ratio of y to x for vapor-liquid equilib calcs;

Parameter
 Q(GnrlE)                Used as a calculation intermediate for initialization only;

*----------------------------------------------------*
*        Equations for General Thermo Equipment
*----------------------------------------------------*

Equations
EqCompMolBalGnrlE(GnrlE,Comp)            Component mole balanace
EqEnrgBalGnrlE(GnrlE)                    Energy balance
EqSumFrac2(ThrmE)                        Sum mole fraction
EqSumFrac3(ThrmE)                        Sum mole fraction
EqVLEThrmE(ThrmE,Comp,Str,Str2)          vapor-liquid equilibrium
EqTRelThrmE(ThrmE,Str,Str2)              thermal equilibrium
EqPRelThrmE(ThrmE,Str,Str2)              pressure equilibrium;

EqCompMolBalGnrlE(GnrlE,Comp)$(NOT InactiveGnrlE(GnrlE))..
  Sum(Str$InGnrlE(GnrlE,Str), Fc(Str,Comp))
         =E= Sum(Str2$OutVGnrlE(GnrlE,Str2), Fc(Str2,Comp))
         +  Sum(Str3$OutLGnrlE(GnrlE,Str3), Fc(Str3,Comp))  ;

EqEnrgBalGnrlE(GnrlE)$(CalcGnrlE(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
  SUM(Str$InGnrlE(GnrlE,Str), F(Str)*H(Str)) + Qin(GnrlE)
         =e= SUM(Str2$OutVGnrlE(GnrlE,Str2), F(Str2)*H(Str2))
         + SUM(Str3$OutLGnrlE(GnrlE, Str3), F(Str3)*H(Str3))
         + Qout(GnrlE);

EqSumFrac2(ThrmE)$(NOT InactiveGnrlE(ThrmE))..
  SUM(Comp, SUM(Str2$OutVGnrlE(ThrmE,Str2), Xc(Str2,Comp))) =e=
  SUM(Comp, SUM(Str3$OutLGnrlE(ThrmE,Str3), Xc(Str3,Comp)))  ;

EqSumFrac3(ThrmE)$(NOT InactiveGnrlE(ThrmE))..
  SUM(Comp, SUM(Str2$OutVGnrlE(ThrmE,Str2), Xc(Str2,Comp))) +
  SUM(Comp, SUM(Str3$OutLGnrlE(ThrmE,Str3), Xc(Str3,Comp))) =e= 2 ;

EqVLEThrmE(ThrmE,Comp,Str,Str2)$(OutVGnrlE(ThrmE,Str) AND OutLGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  Xc(Str,Comp) =e= beta(ThrmE)*K(ThrmE,Comp)*Xc(Str2,Comp);

EqTRelThrmE(ThrmE,Str,Str2)$(OutVGnrlE(ThrmE,Str) AND OutLGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  Tscaled(Str) =e= Tscaled(Str2);

EqPRelThrmE(ThrmE,Str,Str2)$(OutVGnrlE(ThrmE,Str) AND OutLGnrlE(ThrmE,Str2) AND NOT InactiveGnrlE(ThrmE))..
  P(Str) =e= P(Str2);

Model EqThrmE /EqCompMolBalGnrlE, EqEnrgBalGnrlE, EqVLEThrmE, EqTRelThrmE, EqPRelThrmE/;
*EqSumFrac2, EqSumFrac3
