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
    CmprPwr(FeedStr)*(Fc('S15','O2') + Fc('S17','O2')) =e= F(FeedStr)*Ncmpr*gamma/(gamma - 1)*Rsi*300/(1000*115.2)*( P(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1);

Parameter
  QwWeight                               /1/;

variable
  o2pure         Oxygen product purity (mol fraction)
  o2rec          Oxygen recovery for process (mol fraction)
  Z9             Dummy variable for obj func ;

Equations
  eqo2pure                       Oxygen purity calcs
  eqo2rec                        Oxygen recovery calcs
  EqFeedSplit
  EqFeedXc(FeedStr,Comp)         Fix feed composition to match specifications
  obj9
  obj9b;

eqo2pure..
         o2pure*sum(Comp, Fc('S15',Comp) + Fc('S17',Comp)) =e= (Fc('S15','O2')+Fc('S17','O2'));

eqo2rec..
         o2rec*SUM(FeedStr, Fc(FeedStr,'O2') ) =e= (Fc('S15','O2')+Fc('S17','O2'));

obj9..
         Z9  =e= ComplPen*vapPen + ComplPen*liqPen + SUM(FeedStr, CmprPwr(FeedStr))
                 + QwWeight*Qw + Qs;

obj9b..
         Z9 =e= ComplPen*vapPen + ComplPen*liqPen + SUM(FeedStr, CmprPwr(FeedStr))
                 + 5*Qw  - o2rec + Qs;
;

EqFeedSplit(Comp)..
  SUM(FeedStr, Fc(FeedStr,Comp)) =e= 2*FeedFlow(Comp);

EqFeedXc(FeedStr,Comp)$(ord(FeedStr) < card(FeedStr) AND ord(Comp) < card(Comp))..
  Xc(FeedStr,Comp) =e= FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));

Model PurityRecoveryEquations /EqO2Pure, EqO2Rec/;

Model EqObj /PurityRecoveryEquations, EqLiqPen, EqVapPen, EqFeedSplit, EqFeedXc, EqCmprPwr/;

* Section 1 (Simple Thermo, Shortcut Cascade)

* Switched from obj9b to obj9 on May 19th, 2014
Model ObjSec1 /obj9b/;

* Section 3 (CEOS Thermo, Shorcut Cascade)
Model ObjSec3 /obj9/;

* Section 4 (CEOS Thermo, MESH Cascade)
Model ObjSec4 /obj9/;

* Section 6 (CEOS Thermo, MESH Cascade, Heat Exchange Subunits)
Model ObjSec6 /obj9/;
