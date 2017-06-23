***********************
***** ObjFncs.gms *****
***********************

*----------------------------------------------*
**  Variables & Equation for objective function
*----------------------------------------------*

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

Parameter
  QwWeight                               /1/;

variable
  Z9             Dummy variable for obj func ;

Equations
  obj
  DummyEqn(Str);

obj..
          Z9 =e= ComplPen*vapPen + ComplPen*liqPen + Qs + QwWeight*Qw + 0.1*SUM(FeedStr, CmprPwr(FeedStr)) + Tscaled('S4');

DummyEqn(Str)$(1 < 0)..
  F(Str) =g= -1;


Model PurityRecoveryEquations / DummyEqn /;

Model EqObj /EqLiqPen, EqVapPen, EqCmprPwr /;

* Section 1 (Simple Thermo, Shortcut Cascade)
Model ObjSec1 /obj/;

* Section 3 (CEOS Thermo, Shorcut Cascade)
Model ObjSec3 /obj/;

* Section 4 (CEOS Thermo, MESH Cascade)
Model ObjSec4 /obj/;

* Section 6 (CEOS Thermo, MESH Cascade, Heat Exchange Subunits)
Model ObjSec6 /obj/;
