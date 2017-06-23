**********************************
***** HeatExchangerModel.gms *****
**********************************

* This file contains the heat exchanger model. Pinch analysis for heat
* integration is modeled in a different file.

* Heat exchangers are modeled as ThrmE units. As a result, they inherit
* equations from ThrmEModel.gms.

* There are two classifications of heat exchangers:
*   - Heating HX: Units that receive heating
*   - Cooling HX: Unit that receive cooling

* Note: Heating HXs are assosiated with COLD streams for heat integration;
* cooling HXs are assosiated with HOT streams.

* Large heat exchangers are approximated as multiple units. This helps combat
* nonlinearities in heat capacity assosiated with phase changes.

*----------------------------------------*
**  Heat exchanger equations
*----------------------------------------*

Alias (HtEx,HtEx1);

Equations
  EqPresHtEx(HtEx,Str,Str2)              No pressure drop
  EqCoolHtEx(CoolHtEx,Str,Str2)          Temperature decreasing in Cooling HX
  EqHeatHtEx(HeatHtEx,Str,Str2)          Temperature increases in Heating HX
  EqDecreaseLiq(HeatHtEx,Str2)
  EqDecreaseVap(CoolHtEx,Str2)
*  EqSpreadCool(HtExGrps,HtEx)            Distribute heating load equally between paired cooling HX (not used currently)
*  EqSpreadHeat(HtExGrps,HtEx)            Distribute heating load equally between paired heat HX (not used currently);
;

EqPresHtEx(HtEx,Str,Str2)$(InThrmEPres(HtEx,Str) AND OutVHtEx(HtEx,Str2) AND NOT InactiveGnrlE(HtEx))..
  P(Str) =e= P(Str2);

EqCoolHtEx(CoolHtEx,Str,Str2)$(InThrmEPres(CoolHtEx,Str) AND OutVHtEx(CoolHtEx,Str2) AND NOT InactiveGnrlE(CoolHtEx))..
  Tscaled(Str) =g= Tscaled(Str2);

EqHeatHtEx(HeatHtEx,Str,Str2)$(InThrmEPres(HeatHtEx,Str) AND OutVHtEx(HeatHtEx,Str2) AND NOT InactiveGnrlE(HeatHtEx))..
  Tscaled(Str) =l= Tscaled(Str2);

EqDecreaseLiq(HeatHtEx,Str2)$(OutLHtEx(HeatHtEx,Str2) AND NOT InactiveGnrlE(HeatHtEx))..
  SUM(Str$(InHtEx(HeatHtEx,Str) AND (LiqStr(Str) OR FlashLiq(Str))), F(Str)) =g= F(Str2);

EqDecreaseVap(CoolHtEx,Str2)$(OutVHtEx(CoolHtEx,Str2) AND NOT InactiveGnrlE(CoolHtEx))..
  SUM(Str$(InHtEx(CoolHtEx,Str) AND (VapStr(Str) OR FlashVap(Str))), F(Str)) =g= F(Str2);

*EqSpreadCool(HtExGrps,HtEx)$(HtExBig(HtExGrps,HtEx) AND CoolHtEx(HtEx) AND NOT InactiveGnrlE(HtEx))..
*  Qout(HtEx) =e= SUM(HtEx1$HtExBig(HtExGrps,HtEx1), Qout(HtEx1))/2;

*EqSpreadHeat(HtExGrps,HtEx)$(HtExBig(HtExGrps,HtEx) AND HeatHtEx(HtEx) AND NOT InactiveGnrlE(HtEx))..
*  Qin(HtEx) =e= SUM(HtEx1$HtExBig(HtExGrps,HtEx1), Qin(HtEx1))/2;

***** Restrict Qin for cooling HXs and Qout for heating HXs to zero *****
Qin.fx(CoolHtEx) = 0;
Qout.fx(HeatHtEx) = 0;

Model EqHtEx /EqPresHtEx, EqCoolHtEx, EqHeatHtEx, EqDecreaseLiq, EqDecreaseVap/;
* EqCoolHtEx, EqHeatHtEx
* EqSpreadCool, EqSpreadHeat
