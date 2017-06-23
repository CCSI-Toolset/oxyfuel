***********************
***** ObjFncs.gms *****
***********************

* This file contains equations required for the objective function. There are
* several objective functions available in this file.

***** Complementarity Penalties *****

Positive Variables
  liqPen
  vapPen;

Equations
  EqLiqPen               Calculate liquid penalty - used for complementarities
  EqVapPen               Calculate vapor penalty - used for complementarities;

EqVapPen..
 vapPen =e= sum(Str$(FlashVap(Str) AND NOT InactiveStr(Str)), sV(Str)*F(Str));

EqLiqPen..
 liqPen =e= sum(Str$(FlashLiq(Str) AND NOT InactiveStr(Str)), sL(Str)*F(Str));

***** Air Compressor Power Calculations *****

Parameters
  Ncmpr /3/
  gamma /1.4/;

Positive Variable
  CmprPwr(FeedStr);

Equations
  EqCmprPwr(FeedStr);

EqCmprPwr(FeedStr)..
    CmprPwr(FeedStr) =e= Ncmpr*gamma/(gamma - 1)*Rsi*300/(1000*115.2)*( P(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1);

*----------------------------------------------*
**  Variables & Equation for objective function
*----------------------------------------------*

variable
  Z Dummy variable for obj func;

Equations
  obj;

obj..
          Z =e= ComplPen*vapPen + ComplPen*liqPen + Qs + Qw + 0.1*SUM(FeedStr, CmprPwr(FeedStr)) + Tscaled('S4');

Model EqObj /EqLiqPen, EqVapPen, EqCmprPwr, obj /;
