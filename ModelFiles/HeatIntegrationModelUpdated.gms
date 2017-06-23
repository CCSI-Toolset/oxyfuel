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

*Positive Variables
*  Tin(GnrlE)                      Inlet temperature
*  Tout(GnrlE)                     Outlet temperature

*  TPcnd(StrBase)                 Pinch candidate temperatures
*  Qs                             Steam (hot) utility load
*  Qw                             Water (cold) utility load;

*Equations
*  EqQw                           Calculate cooling water demand from energy balance + Qs
*  EqMaxDT_cool(HtEx)             Specify max change in temperature for cooling HXs
*  EqMaxDT_heat(HtEx)             Specify max change in temperature for heating HXs
*  EqMaxDT_cond(TCond)            Specify max change in temperature for condensers
*  EqMaxDT_reb(PReb)              Specify max change in temperature for reboilers
*  ;

*EqQw..
*  Qw =e= Qs + SUM(GnrlE$(NOT InactiveGnrlE(GnrlE) AND HeatInteg(GnrlE)), Qout(GnrlE) - Qin(GnrlE));

*EqMaxDT_cool(HtEx)$CoolHtEx(HtEx)..
*  SUM(Str$InHtExOne(HtEx,Str), T(Str)) - SUM(Str2$OutVHtEx(HtEx,Str2), T(Str2)) =l= MaxDT;

*EqMaxDT_heat(HtEx)$HeatHtEx(HtEx)..
*  SUM(Str2$OutVHtEx(HtEx,Str2), T(Str2)) - SUM(Str$InHtExOne(HtEx,Str), T(Str)) =l= MaxDT;

*EqMaxDT_cond(TCond)..
*  SUM(Str$InTCond(TCond,Str), T(Str)) - SUM(Str2$OutTCond(TCond,Str2), T(Str2)) =l= MaxDT;

*EqMaxDT_reb(PReb)..
*  SUM(Str2$OutVPReb(PReb,Str2), T(Str2)) - SUM(Str$InPReb(PReb,Str), T(Str)) =l= MaxDT;

*----------------------------------------*
**  Heat integration equations
*----------------------------------------*

Alias (HtEx, HtEx2), (PReb, PReb2), (Tcond, Tcond2);

Variables
  QAhU(Str,HIZone)                   Available heating above the pinch candidate (in set Str)
  QAcU(Str,HIZone)                   Availabe cooling above the pinch candidate (in set Str);

Positive Variables
  FCp(GnrlE)                     Flow times heat capacity;

Equations
  EqTin_heat(GnrlE,Str,Str2)      Calculate Tin for heating units
  EqTout_heat(GnrlE,Str,Str2)     Calculate Tout for heating units
  EqTin_cool(GnrlE,Str,Str2)      Calculate Tin for cooling units
  EqTout_cool(GnrlE,Str,Str2)     Calculate Tout for cooling units


  EqCalcFCp_heat(GnrlE)           Calculate FCp for heating units
  EqCalcFCp_cool(GnrlE)           Calculate FCp for cooling units

  EqTPcndHotU(GnrlE,Str,HIZone)       Assign pinch candidate temperatures for cooling units (inlet)
  EqTPcndColdU(GnrlE,Str,HIZone)      Assign pinch candidate temperatures for heating units (inlet)

  EqTPcndHot2U(GnrlE,Str,HIZone)      Assign pinch candidate temperatures for cooling units (outlet)
  EqTPcndCold2U(GnrlE,Str,HIZone)     Assign pinch candidate temperatures for heating units (outlet)

  EqQAhU(Str,HIZone)                 Calculate heating above pinch for all candidates
  EqQAcU(Str,HIZone)                 Calculate cooling above pinch for all candidates

  EqQAhUTrick(Str,HIZone)            Calculate heating above pinch for all candidates
  EqQAcUTrick(Str,HIZone)            Calculate cooling above pinch for all candidates

  EqQsU(Str,HIZone)                  Calculate steam demand for all candidates from QAc and QAh;

***** Cold Streams *****
EqTin_heat(GnrlE,Str,Str2)$(InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
* Original
  Tin(GnrlE) =e= Tref*Tscaled(Str);

* Alternate Smooth
*  Tin(GnrlE) =e= T(Str2) - smmax(T(Str2) - T(Str) + alpha) + alpha;

* Constant Bump
*  Tin(GnrlE) =e= T(Str) + alpha;

EqTout_heat(GnrlE,Str,Str2)$(InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
* Original
  Tout(GnrlE) =e= smmax(Tref*Tscaled(Str2) - Tref*Tscaled(Str) + alpha) + Tref*Tscaled(Str) - alpha;

* Alternate Smooth/Constant Bump
*  Tout(GnrlE) =e= T(Str2);

***** Hot Streams *****
EqTin_cool(GnrlE,Str,Str2)$(InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
* Original
  Tin(GnrlE) =e= smmax(Tref*Tscaled(Str) - Tref*Tscaled(Str2) + alpha) + Tref*Tscaled(Str2) - alpha;

* Constant Bump
*  Tin(GnrlE) =e= T(Str);

* Alternate Smooth
*  Tin(GnrlE) =e= T(Str2) + smmax(T(Str) - T(Str2) + alpha) - alpha;

EqTout_cool(GnrlE,Str,Str2)$(InOneGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2) AND CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
* Constant Bump
*  Tout(GnrlE) =e= T(Str2) - alpha;

* Original & Alternate Smooth
  Tout(GnrlE) =e= Tref*Tscaled(Str2);

EqCalcFCp_heat(GnrlE)$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
  FCp(GnrlE)*(Tout(GnrlE) - Tin(GnrlE)) =e= Qin(GnrlE);

EqCalcFCp_cool(GnrlE)$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE))..
  FCp(GnrlE)*(Tin(GnrlE) - Tout(GnrlE)) =e= Qout(GnrlE);


EqTPcndHotU(GnrlE,Str,HIZone)$(CEqp(GnrlE) AND InOneGnrlE(GnrlE,Str) AND NOT InactiveGnrlE(GnrlE) AND HIMap(GnrlE, HIZone))..
  TPcnd(Str,HIZone) =e= Tin(GnrlE);

EqTPcndHot2U(GnrlE,Str,HIZone)$(CEqp(GnrlE) AND OutLGnrlE(GnrlE,Str) AND OutPStr(Str,HIZone) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(Str) AND HIMap(GnrlE, HIZone))..
  TPcnd(Str,HIZone) =e= Tout(GnrlE);

EqTPcndColdU(GnrlE,Str,HIZone)$(HEqp(GnrlE) AND InOneGnrlE(GnrlE,Str) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(Str) AND HIMap(GnrlE, HIZone))..
  TPcnd(Str,HIZone) =e= Tin(GnrlE) + HRAT(HIZone);

EqTPcndCold2U(GnrlE, Str,HIZone)$(HEqp(GnrlE) AND OutLGnrlE(GnrlE,Str) AND OutPStr(Str,HIZone) AND NOT InactiveGnrlE(GnrlE) AND NOT InactiveStr(Str) AND HIMap(GnrlE, HIZone))..
  TPcnd(Str,HIZone) =e= Tout(GnrlE) + HRAT(HIZone);

EqQAhU(Str,HIZone)$(PStr(Str,HIZone) AND NOT InactiveStr(Str))..
  QAhU(Str,HIZone) =e= SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND HIMap(GnrlE, HIZone)), FCp(GnrlE)*Tref*smmax(Tin(GnrlE) - TPcnd(Str,HIZone)) - FCp(GnrlE)*Tref*smmax( Tout(GnrlE) - TPcnd(Str,HIZone) ) );

EqQAhUTrick(Str,HIZone)$(PStr(Str,HIZone) AND NOT InactiveStr(Str))..
  QAhU(Str,HIZone) =e= SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,Str) AND NOT OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)),
                         FCp(GnrlE)*smmax(Tin(GnrlE) - TPcnd(Str,HIZone)) - FCp(GnrlE)*smmax( Tout(GnrlE) - TPcnd(Str,HIZone) ) )
                   + SUM(GnrlE$(CEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)), Qout(GnrlE));

EqQAcU(Str,HIZone)$(PStr(Str,HIZone) AND NOT InactiveStr(Str))..
  QAcU(Str,HIZone) =e= SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND HIMap(GnrlE, HIZone)), FCp(GnrlE)*smmax( Tout(GnrlE) - TPcnd(Str,HIZone) + HRAT(HIZone)) - FCp(GnrlE)*Tref*smmax(Tin(GnrlE) - TPcnd(Str,HIZone) + HRAT(HIZone) ) );

EqQAcUTrick(Str,HIZone)$(PStr(Str,HIZone) AND NOT InactiveStr(Str))..
  QAcU(Str,HIZone) =e= SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND NOT InOneGnrlE(GnrlE,Str) AND NOT OutLGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)),
                         FCp(GnrlE)*smmax( Tout(GnrlE) - TPcnd(Str,HIZone) + HRAT(HIZone)) - FCp(GnrlE)*smmax(Tin(GnrlE) - TPcnd(Str,HIZone) + HRAT(HIZone) ) )
                   + SUM(GnrlE$(HEqp(GnrlE) AND NOT InactiveGnrlE(GnrlE) AND InOneGnrlE(GnrlE,Str) AND HIMap(GnrlE, HIZone)), Qin(GnrlE));



EqQsU(Str,HIZone)$(PStr(Str, HIZone) AND NOT InactiveStr(Str))..
  QsZ(HIZone) =g= QAcU(Str, HIZone) - QAhU(Str, HIZone);

Model HeatIntegUpdated /EqTin_heat, EqTout_heat, EqTin_cool, EqTout_cool,
                         EqCalcFCp_cool, EqCalcFCp_heat,
                         EqTPcndHotU,EqTPcndColdU,
                         EqTPcndHot2U, EqTPcndCold2U,
                         EqQAhUTrick, EqQAcUTrick, EqQsU, EqQwZ, EqQsTotal, EqQwTotal, EqMaxDT_cond, EqMaxDT_reb/;
