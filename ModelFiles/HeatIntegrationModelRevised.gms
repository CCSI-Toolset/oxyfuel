************************************
***** HeatIntegrationModel.gms *****
************************************

* This file contains the heat integration model, which is based on the
* optimization friendly pinch location formulation of Duran & Grossmann.
* In summary the available heating and cooling above every pinch point is
* calculated. They are used to calculate the minimium required heating and
* cooling utilities required while ensure a minimium approach temperature (HRAT)
* is maintained. For additional information, please see the following papers:
*
* "Simultaneous optimization and heat integration of chemical processes" by
* Duran & Grossmann (1986)
*
* "Modeling multistream heat exchangers with and without phase changes for
* simultaneous optimization and heat integration" by Kamath, Biegler
* & Grossmann (2012)
*
* In this model all heat exchangers, condensers and reboilers in the
* flowsheet are considered for heat integration.
*
* In this model smoothed max operators are used both for available duty
* calculates and to "bump" temperatures with a phase change. The "bump" ensures
* there is a temperature difference of at least alpha for heat capacity
* calculations.
*
* Here is the mapping from units to COLD/HOT streams:
*   - Cooling HX --> HOT stream
*   - Condenser --> HOT stream
*   - Heating HX --> COLD stream
*   - Reboiler --> COLD stream

*----------------------------------------*
**  Heat integration equations
*----------------------------------------*

Variables
  QAhR(StrBase,GnrlE)               Available heating above the pinch candidate
  QAcR(StrBase,GnrlE)               Availabe cooling above the pinch candidate
  a(StrBase,GnrlE)                a for heat integration equipment
;

a.lo(StrBase,HeatInteg) = 0;
a.up(StrBase,HeatInteg) = 1;


Positive Variables
  TPcnd(StrBase)                 Pinch candidate temperatures
  Qs                             Steam (hot) utility load
  Qw                             Water (cold) utility load;

Equations
  EqTPcndHotR(GnrlE,StrBase)       Assign pinch candidate temperatures for cooling equipment (inlet)
  EqTPcndColdR(GnrlE,StrBase)      Assign pinch candidate temperatures for heating equipment (inlet)
  EqTPcndHot2R(GnrlE,StrBase)      Assign pinch candidate temperatures for cooling equipment (outlet)
  EqTPcndCold2R(GnrlE,StrBase)     Assign pinch candidate temperatures for heating equipment (outlet)

  EqQAhR(StrBase,GnrlE,Str,Str2)       Calculate heating above pinch for all candidates
  EqQhRExtra(StrBase,GnrlE)

  EqQAcR(StrBase,GnrlE,Str,Str2)       Calculate cooling above pinch for all candidates

  EqQAhRSpecialIn(StrBase,GnrlE)
  EqQAhRSpecialOut(StrBase,GnrlE)

  EqQAcRSpecialIn(StrBase,GnrlE)
  EqQAcRSpecialOut(StrBase,GnrlE)

  EqQsR(StrBase)                  Calculate steam demand for all candidates from QAc and QAh
  ;

EqTPcndHotR(GnrlE,StrBase)$(CEqp(GnrlE) AND InOneGnrlE(GnrlE,StrBase) AND NOT InactiveGnrlE(GnrlE))..
  TPcnd(StrBase) =e= T(StrBase);

EqTPcndHot2R(GnrlE, StrBase)$( CEqp(GnrlE) AND OutLGnrlE(GnrlE,StrBase) AND OutPStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= T(StrBase) + alpha;

EqTPcndColdR(GnrlE, StrBase)$(HEqp(GnrlE) AND InOneGnrlE(GnrlE,StrBase) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= T(StrBase) + HRAT + alpha;

EqTPcndCold2R(GnrlE, StrBase)$(HEqp(GnrlE) AND OutLGnrlE(GnrlE,StrBase) AND OutPStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= T(StrBase) + HRAT;

* Hot Streams
EqQAhR(StrBase,GnrlE,Str,Str2)$(CEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,StrBase) AND NOT OutLGnrlE(GnrlE,StrBase))..
  (T(Str) - T(Str2) - alpha)*QAhR(StrBase,GnrlE) =l= a(StrBase,GnrlE)*Qout(GnrlE)*(T(Str) - TPcnd(StrBase));

EqQhRExtra(StrBase,GnrlE)$(CEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,StrBase) AND NOT OutLGnrlE(GnrlE,StrBase))..
  QAhR(StrBase,GnrlE) =l= Qout(GnrlE);

EqQAhRSpecialIn(StrBase,GnrlE)$(CEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND InOneGnrlE(GnrlE,StrBase))..
  QAhR(StrBase,GnrlE) =e= 0;

EqQAhRSpecialOut(StrBase,GnrlE)$(CEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND OutLGnrlE(GnrlE,StrBase))..
  QAhR(StrBase,GnrlE) =e= Qout(GnrlE);

* Cold Stream
EqQAcR(StrBase,GnrlE,Str,Str2)$(HEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,StrBase) AND NOT OutLGnrlE(GnrlE,StrBase))..
  (T(Str2) - T(Str)  - alpha)*QAcR(StrBase,GnrlE) =g= Qin(GnrlE)*a(StrBase,GnrlE)*(T(Str2) - TPcnd(StrBase) + HRAT) + (1-a(StrBase,GnrlE))*Qin(GnrlE)*(T(Str2) - T(Str) - alpha);

EqQAcRSpecialIn(StrBase,GnrlE)$(HEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND InOneGnrlE(GnrlE,StrBase))..
  QAcR(StrBase,GnrlE) =e= Qin(GnrlE);

EqQAcRSpecialOut(StrBase,GnrlE)$(HEqp(GnrlE) AND PStr(StrBase) AND NOT InactiveStr(StrBase) AND NOT InactiveGnrlE(GnrlE) AND OutLGnrlE(GnrlE,StrBase))..
  QAcR(StrBase,GnrlE) =e= 0;


QAcR.lo(StrBase,HEqp) = 0;

QAcR.up(StrBase,HEqp) = 100;

QAhR.lo(StrBase,CEqp) = -10;
QAhR.up(StrBase,CEqp) = 1E4;

Qs.up = 1E5;
Qw.up = 1E5;


EqQsR(StrBase)$(PStr(StrBase) AND NOT InactiveStr(StrBase))..
  Qs =g= SUM(HEqp$(NOT InactiveGnrlE(HEqp)), QAcR(StrBase, HEqp))
       - SUM(CEqp$(NOT InactiveGnrlE(CEqp)), QAhR(StrBase, CEqp));

Model HeatIntegRevised /EqTPcndHotR,EqTPcndColdR,
                         EqTPcndHot2R, EqTPcndCold2R,
                         EqQAhR, EqQAcR,
                         EqQAhRSpecialIn, EqQAhRSpecialOut,
                         EqQAcRSpecialIn, EqQAcRSpecialOut,
                         EqQhRExtra,
                         EqQsR, EqQw,
                         EqMaxDT_cond, EqMaxDT_reb/;

* Parameter for initialization
Scalar Asmall            /0/;
Scalar Alarge            /1/;
