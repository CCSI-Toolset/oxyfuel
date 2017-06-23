***********************
***** ObjFncs.gms *****
***********************

* This file contains equations required for the objective function. There are
* several objective functions available in this file.

***** Complementarity Penalties *****

Positive Variables
  liqPen
  vapPen
  heatPen
;

Equations
  EqLiqPen               Calculate liquid penalty - used for complementarities
  EqVapPen               Calculate vapor penalty - used for complementarities
  EqHeatPen
;

EqVapPen..
 vapPen =e= sum(Str$(FlashVap(Str) AND NOT InactiveStr(Str)), sV(Str)*F(Str));

EqLiqPen..
 liqPen =e= sum(Str$(FlashLiq(Str) AND NOT InactiveStr(Str)), sL(Str)*F(Str));

EqHeatPen..
  heatPen =e= 0;

*----------------------------------------------*
**  Variables & Equation for objective function
*----------------------------------------------*

Parameter
  QwWeight                               /1/;

variable
  o2pure         Oxygen product purity (mol fraction)
  o2rec          Oxygen recovery for process (mol fraction)
  Z9             Dummy variable for obj func ;
Equations
  eqo2pure                       Oxygen purity calcs
  eqo2rec                        Oxygen recovery calcs
  EqFeedXc(FeedStr,Comp)         Fix feed composition to match specifications
  obj9
  obj9b;

eqo2pure..
         o2pure =e= Xc('S10','O2');

eqo2rec..
         o2rec*SUM(FeedStr, Fc(FeedStr,'O2') ) =e= Fc('S10','O2');

obj9..
         Z9  =e= ComplPen*vapPen + ComplPen*liqPen
                 + QwWeight*Qw + Qs + heatPen;

obj9b..
         Z9 =e= ComplPen*vapPen + ComplPen*liqPen
                 + 5*Qw  - o2rec  + heatPen;
;

EqFeedXc(FeedStr,Comp)$(ord(FeedStr) < card(FeedStr) AND ord(Comp) < card(Comp))..
  Xc(FeedStr,Comp) =e= FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));

Model EqObj /EqO2Pure, EqO2Rec, EqLiqPen, EqVapPen, EqFeedXc, EqHeatPen/;
