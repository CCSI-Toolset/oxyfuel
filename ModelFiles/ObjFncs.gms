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
  heatPen =e= Qin('HX15')*Qout('HX16');

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

*----------------------------------------------*
**  Variables & Equation for objective function
*----------------------------------------------*

Parameter
  QwWeight                               /1/;

variable
  o2pure         Oxygen product purity (mol fraction)
  o2rec          Oxygen recovery for process (mol fraction)
  Z2             Dummy variable for obj func
  Z3             Dummy variable for obj func
  Z4             Dummy variable for obj func
  Z5             Dummy variable for obj func
  Z7             Dummy variable for obj func
  Z8             Dummy variable for obj func
  Z9             Dummy variable for obj func
  Z10            Dummy variable for obj func;

Equations
  eqo2pure                       Oxygen purity calcs
  eqo2rec                        Oxygen recovery calcs
  EqFeedSplit(Comp)              Ensure feed sums to 2 mole per time
  EqFeedXc(FeedStr,Comp)         Fix feed composition to match specifications
  obj
  obj2
  obj3
  obj4
  obj5
  obj7
  obj8
  obj9
  obj9b
  obj9c
  objHX1
  objHX2
  obj10;

eqo2pure..
         o2pure*sum(Comp, Fc('S15',Comp) + Fc('S17',Comp)) =e= (Fc('S15','O2')+Fc('S17','O2'));

eqo2rec..
         o2rec*SUM(FeedStr, Fc(FeedStr,'O2') ) =e= (Fc('S15','O2')+Fc('S17','O2'));


* Replaced 0.9 with o2pureSpec
obj..
*         Z =e= (o2pure - o2pureSpec)*(o2pure - o2pureSpec) + vapPen + liqPen - o2rec + Qw + Qs + 0.1*SUM(Valve, Qin(Valve) + Qout(Valve));
          Z =e= ComplPen*vapPen + ComplPen*liqPen + 0.1*Qs + 0.1*Qw;
*         Z  =e= ComplPen*vapPen + ComplPen*liqPen
*                 + Qw - 0.1*Qs + SUM(FeedStr, CmprPwr(FeedStr));

objHX1..
*        Z =e= ComplPen*vapPen + ComplPen*liqPen + Qs + Qw;
         Z  =e= ComplPen*vapPen + ComplPen*liqPen
                 + Qw - 0.1*Qs + SUM(FeedStr, CmprPwr(FeedStr));

objHX2..
        Z =e= ComplPen*vapPen + ComplPen*liqPen - o2rec + 0.1*Qs + Qw;
*         Z =e= vapPen + liqPen - o2rec + Qs + Qw;

obj2..
*         Z2 =e= ComplPen*vapPen + ComplPen*liqPen
*                 + sum(ThrmE, Qout(ThrmE) - 0.5*Qin(ThrmE)) - 0.5*sum(HtEx, Qin(HtEx))
*                + sum(TCond, QTCond(TCond)) + 0.01*SUM(FeedStr, P(FeedStr));

          Z2 =e= ComplPen*vapPen + ComplPen*liqPen + Qw + Qs - o2rec;

*         Z2 =e= ComplPen*vapPen + ComplPen*liqPen + Qw + Qs
*                 + 0.1*SUM(FeedStr, P(FeedStr));


obj3..
         Z3 =e= ComplPen*vapPen + ComplPen*liqPen -o2pure;

obj4..
         Z4 =e= ComplPen*vapPen + ComplPen*liqPen - o2pure + Qs + Qw;

obj5..   Z5 =e= ComplPen*vapPen + ComplPen*liqPen
                 + sum(GnrlE, Qin(GnrlE) + Qout(GnrlE))
                 - F('SF0')*T('SF0') - F('SF10')*T('SF10');

obj7..
         Z7 =e= ComplPen*vapPen + ComplPen*liqPen + Qs + Qw;

obj8..
         Z8*o2rec =e= ComplPen*vapPen + ComplPen*liqPen + Qs + Qw;

obj9..
         Z9  =e= ComplPen*vapPen + ComplPen*liqPen
                 + QwWeight*Qw + Qs + SUM(FeedStr, CmprPwr(FeedStr)) + heatPen;

obj9b..
         Z9 =e= ComplPen*vapPen + ComplPen*liqPen
                 + 5*Qw + SUM(FeedStr, CmprPwr(FeedStr))  - o2rec  + heatPen;

obj9c..
         Z9 =e= ComplPen*vapPen + ComplPen*liqPen
                 + 5*Qw + 0.1*Qs + SUM(FeedStr, CmprPwr(FeedStr)) + o2pure;

obj10..   Z10 =e= ComplPen*vapPen + ComplPen*liqPen
                 + sum(GnrlE, Qin(GnrlE) + Qout(GnrlE));

EqFeedSplit(Comp)..
  SUM(FeedStr, Fc(FeedStr,Comp)) =e= 2*FeedFlow(Comp);

EqFeedXc(FeedStr,Comp)$(ord(FeedStr) < card(FeedStr) AND ord(Comp) < card(Comp))..
  Xc(FeedStr,Comp) =e= FeedFlow(Comp)/SUM(Comp2, FeedFlow(Comp2));

Model EqObj /EqO2Pure, EqO2Rec, EqLiqPen, EqVapPen, EqFeedSplit, EqFeedXc, EqHeatPen/;
