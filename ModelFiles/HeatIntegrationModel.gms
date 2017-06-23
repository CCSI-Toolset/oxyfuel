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

Alias (HtEx, HtEx2), (PReb, PReb2), (Tcond, Tcond2), (PStr, PStr2);

Variables
  QAh(StrBase)                   Available heating above the pinch candidate (in set StrBase)
  QAc(StrBase)                   Availabe cooling above the pinch candidate (in set StrBase);

Positive Variables
  Tin(HtEx)                      Inlet temperature for HtEx
  Tout(HtEx)                     Outlet temperature for HtEx
  Tin_r(PReb)                    Inlet temperature for PReb
  Tout_r(PReb)                   Outlet temperature for PReb
  Tin_c(TCond)                   Inlet temperature for TCond
  Tout_c(TCond)                  Outlet temperature for TCond
  FCp(HtEx)                      Flow times heat capacity for HtEx
  FCp_cond(Tcond)                Flow times heat capacity for TCond
  FCp_reb(Preb)                  Flow times heat capacity for PReb
  TPcnd(StrBase)                 Pinch candidate temperatures
  Qs                             Steam (hot) utility load
  Qw                             Water (cold) utility load;

Equations
  EqTin_heat(HtEx,Str,Str2)      Calculate Tin for heating HXs
  EqTout_heat(HtEx,Str,Str2)     Calculate Tout for heating HXs
  EqTin_cool(HtEx,Str,Str2)      Calculate Tin for cooling HXs
  EqTout_cool(HtEx,Str,Str2)     Calculate Tout for cooling HXs
  EqTin_r(PReb,Str,Str2)         Calculate Tin for partial reboilers
  EqTout_r(PReb,Str,Str2)        Calculate Tout for partial reboilers
  EqTin_c(TCond,Str,Str2)        Calculate Tin for total condensers
  EqTout_c(TCond,Str,Str2)       Calculate Tout for total condensers
  EqCalcFCp_heat(HtEx)           Calculate FCp for heating HXs
  EqCalcFCp_cool(HtEx)           Calculate FCp for cooling HXs
  EqCalcFCp_cond(Tcond)          Calculate FCp for total condensers
  EqCalcFCp_reb(PReb)            Calculate FCp for partial reboilers
  EqTPcndHot(HtEx,StrBase)       Assign pinch candidate temperatures for cooling HXs (inlet)
  EqTPcnd_cond(Tcond,StrBase)    Assign pinch candidate temperatures for total condensers (inlet)
  EqTPcndCold(HtEx,StrBase)      Assign pinch candidate temperatures for heating HXs (inlet)
  EqTPcnd_reb(PReb,StrBase)      Assign pinch candidate temperatures for partial reboilers (inlet)
  EqTPcndHot2(HtEx,StrBase)      Assign pinch candidate temperatures for cooling HXs (outlet)
  EqTPcnd_cond2(Tcond,StrBase)   Assign pinch candidate temperatures for total condensers (outlet)
  EqTPcndCold2(HtEx,StrBase)     Assign pinch candidate temperatures for heating HXs (outlet)
  EqTPcnd_reb2(PReb,StrBase)     Assign pinch candidate temperatures for partial reboilers (outlet)
  EqQAh(StrBase)                 Calculate heating above pinch for all candidates
  EqQAc(StrBase)                 Calculate cooling above pinch for all candidates
  EqQs(StrBase)                  Calculate steam demand for all candidates from QAc and QAh
  EqQw                           Calculate cooling water demand from energy balance + Qs
  EqMaxDT_cool(HtEx)             Specify max change in temperature for cooling HXs
  EqMaxDT_heat(HtEx)             Specify max change in temperature for heating HXs
  EqMaxDT_cond(TCond)            Specify max change in temperature for condensers
  EqMaxDT_reb(PReb)              Specify max change in temperature for reboilers;

EqTin_heat(HtEx,Str,Str2)$(InHtExOne(HtEx,Str) AND OutVHtEx(HtEx,Str2) AND HeatHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  Tin(HtEx) =e= T(Str);

EqTout_heat(HtEx,Str,Str2)$(InHtExOne(HtEx,Str) AND OutVHtEx(HtEx,Str2) AND HeatHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  Tout(HtEx) =e= smmax(T(Str2) - T(Str) + alpha) + T(Str) - alpha;

EqTin_cool(HtEx,Str,Str2)$(InHtExOne(HtEx,Str) AND OutVHtEx(HtEx,Str2) AND CoolHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  Tin(HtEx) =e= smmax(T(Str) - T(Str2) + alpha) + T(Str2) - alpha;

EqTout_cool(HtEx,Str,Str2)$(InHtExOne(HtEx,Str) AND OutVHtEx(HtEx,Str2) AND CoolHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  Tout(HtEx) =e= T(Str2);

EqTin_r(PReb,Str,Str2)$(InPReb(PReb,Str) AND OutVPReb(PReb,Str2) AND NOT InactiveThrmE(PReb))..
  Tin_r(PReb) =e= T(Str);

EqTout_r(PReb,Str,Str2)$(InPReb(PReb,Str) AND OutVPReb(PReb,Str2) AND NOT InactiveThrmE(PReb))..
  Tout_r(PReb) =e= smmax(T(Str2) - T(Str) + alpha) + T(Str) - alpha;

EqTin_c(TCond,Str,Str2)$(InTCond(TCond,Str) AND OutTCond(TCond,Str2))..
  Tin_c(TCond) =e= smmax(T(Str) - T(Str2) + alpha) + T(Str2) - alpha;

EqTout_c(TCond,Str,Str2)$(InTCond(TCond,Str) AND OutTCond(TCond,Str2))..
  Tout_c(TCond) =e= T(Str2);

EqCalcFCp_heat(HtEx)$(HeatHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  FCp(HtEx)*(Tout(HtEx) - Tin(HtEx)) =e= Qin(HtEx);

EqCalcFCp_cool(HtEx)$(CoolHtEx(HtEx) AND NOT InactiveThrmE(HtEx))..
  FCp(HtEx)*(Tin(HtEx) - Tout(HtEx)) =e= Qout(HtEx);

EqCalcFCp_cond(TCond)..
  FCp_cond(TCond)*(Tin_c(TCond) - Tout_c(TCond)) =e= QTcond(Tcond);

EqCalcFCp_reb(PReb)$(NOT InactiveThrmE(PReb))..
  FCp_reb(PReb)*(Tout_r(PReb) - Tin_r(PReb)) =e= Qin(PReb);

EqTPcndHot(HtEx,StrBase)$(CoolHtEx(HtEx) AND InHtExOne(HtEx,StrBase) AND NOT InactiveThrmE(HtEx))..
  TPcnd(StrBase) =e= Tin(HtEx);

EqTPcndHot2(HtEx, StrBase)$( CoolHtEx(HtEx) AND OutVHtEx(HtEx,StrBase) AND OutPStr(StrBase) AND NOT InactiveThrmE(HtEx) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= Tout(HtEx);

EqTPcnd_cond(Tcond,StrBase)$(InTCond(TCond,StrBase))..
  TPcnd(StrBase) =e= Tin_c(Tcond);

EqTPcnd_cond2(Tcond,StrBase)$(OutTCond(TCond,StrBase) AND OutPStr(StrBase) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= Tout_c(Tcond);

EqTPcndCold(HtEx,StrBase)$(HeatHtEx(HtEx) AND InHtExOne(HtEx,StrBase) AND NOT InactiveThrmE(HtEx) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= Tin(HtEx) + HRAT;

EqTPcndCold2(HtEx, StrBase)$(HeatHtEx(HtEx) AND OutVHtEx(HtEx,StrBase) AND OutPStr(StrBase) AND NOT InactiveThrmE(HtEx) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= Tout(HtEx) + HRAT;

EqTPcnd_reb(PReb, StrBase)$InPReb(PReb,StrBase)..
  TPcnd(StrBase) =e= Tin_r(PReb) + HRAT;

EqTPcnd_reb2(PReb, StrBase)$(OutVPReb(PReb,StrBase) AND OutPStr(StrBase) AND NOT InactiveStr(StrBase))..
  TPcnd(StrBase) =e= Tout_r(PReb) + HRAT;

EqQAh(StrBase)$(PStr(StrBase) AND NOT InactiveStr(StrBase))..
  QAh(StrBase) =e= SUM(HtEx$(CoolHtEx(HtEx) AND NOT InactiveThrmE(HtEx)), FCp(HtEx)*smmax(Tin(HtEx) - TPcnd(StrBase)) - FCp(HtEx)*smmax( Tout(HtEx) - TPcnd(StrBase) ) )
            +  SUM(TCond, FCp_cond(Tcond)*smmax(Tin_c(Tcond) - TPcnd(StrBase)) - FCp_cond(Tcond)*smmax( Tout_c(Tcond) - TPcnd(StrBase) ) );

EqQAc(StrBase)$(PStr(StrBase) AND NOT InactiveStr(StrBase))..
  QAc(StrBase) =e= SUM(HtEx$(HeatHtEx(HtEx) AND NOT InactiveThrmE(HtEx)), FCp(HtEx)*smmax( Tout(HtEx) - TPcnd(StrBase) + HRAT) - FCp(HtEx)*smmax(Tin(HtEx) - TPcnd(StrBase) + HRAT ) )
            +  SUM(PReb$(NOT InactiveThrmE(PReb)), FCp_reb(PReb)*smmax( Tout_r(PReb) - TPcnd(StrBase) + HRAT) - FCp_reb(PReb)*smmax(Tin_r(PReb) - TPcnd(StrBase) + HRAT ) );

EqQs(StrBase)$(PStr(StrBase) AND NOT InactiveStr(StrBase))..
  Qs =g= QAc(StrBase) - QAh(StrBase);

EqQw..
  Qw =e= Qs + SUM(HtEx$(NOT InactiveThrmE(HtEx)), Qout(HtEx) - Qin(HtEx))
            + SUM(TCond, QTCond(TCond))
            - SUM(PReb$(NOT InactiveThrmE(PReb)), Qin(PReb));

EqMaxDT_cool(HtEx)$CoolHtEx(HtEx)..
  SUM(Str$InHtExOne(HtEx,Str), T(Str)) - SUM(Str2$OutVHtEx(HtEx,Str2), T(Str2)) =l= MaxDT;

EqMaxDT_heat(HtEx)$HeatHtEx(HtEx)..
  SUM(Str2$OutVHtEx(HtEx,Str2), T(Str2)) - SUM(Str$InHtExOne(HtEx,Str), T(Str)) =l= MaxDT;

EqMaxDT_cond(TCond)..
  SUM(Str$InTCond(TCond,Str), T(Str)) - SUM(Str2$OutTCond(TCond,Str2), T(Str2)) =l= MaxDT;

EqMaxDT_reb(PReb)..
  SUM(Str2$OutVPReb(PReb,Str2), T(Str2)) - SUM(Str$InPReb(PReb,Str), T(Str)) =l= MaxDT;

Model HeatIntegEqns /EqTin_heat, EqTout_heat, EqTin_cool, EqTout_cool,
                         EqCalcFCp_cool, EqCalcFCp_heat,
                         EqTPcndHot,EqTPcndCold,
                         EqTPcndHot2, EqTPcndCold2,
                         EqQAh, EqQAc, EqQs, EqQw/;

Model HeatIntegEqnsExtra /EqTin_r, EqTout_r, EqTin_c, EqTout_c,
                          EqCalcFCp_cond, EqCalcFCp_reb,
                          EqTPcnd_cond, EqTPcnd_reb,
                          EqTPcnd_cond2, EqTPcnd_reb2,
                          EqMaxDT_cond, EqMaxDT_reb/;
