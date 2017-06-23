********************************************
***** Zero Flowrate Degeneracy Example *****
********************************************

* Create by Alex Dowling
* Last modified on May 27th, 2014

$onempty

Sets
  s              streams         /S2*S7/
  c              components      /A*B/
  u              units           /1*2/
  Inlet(u,s)     inlet streams   /(1.S2), (2.S5)/
  OutletV(u,s)   outlet vapor    /(1.S3), (2.S6)/
  OutletL(u,s)   outlet liquid   /(1.S4), (2.S7)/
  FeedStr(s)     feed streams    /S2, S5/
*  DisableMS(u,c) for demo        /(1.A)/
  DisableMS(u,c) for demo        / /
;

Alias (c, c2);

Table
  K(u,c)         equilibrium coefficients
                 A               B
1                1.088           0.9
2                1.099           0.9 ;


Positive Variables
  F(s)           overall stream flowrate (moles per time)
  Fc(s,c)        component molar flowrate (moles per time)
  x(s,c)         mole fraction (vapor and liquid);

Parameters
  FeedFlow(c)    feed molar flowrate     /A = 0.55, B = 0.45/;

***** General Unit Models *****

Equations
  EqUnitComponentMoleBalance(u,c)
  EqUnitMoleBalance(u)
  EqUnitVLE(u,c)
  EqUnitSummation(u)
  EqStrMoleFrac(s,c);

EqUnitComponentMoleBalance(u,c)$(NOT DisableMS(u,c))..
  SUM(s$Inlet(u,s), Fc(s,c)) =e= SUM(s$(OutletV(u,s) OR OutletL(u,s)), Fc(s,c));

EqUnitMoleBalance(u)..
  SUM(s$Inlet(u,s), F(s)) =e= SUM(s$(OutletV(u,s) OR OutletL(u,s)), F(s));

EqUnitVLE(u,c)..
  SUM(s$OutletV(u,s), x(s,c)) =e= K(u,c)*SUM(s$OutletL(u,s), x(s,c));

EqUnitSummation(u)..
  SUM(c, SUM(s$OutletV(u,s), x(s,c)) - SUM(s$OutletL(u,s), x(s,c))) =e= 0;

EqStrMoleFrac(s,c)..
  F(s)*x(s,c) =e= Fc(s,c);

***** Flowsheet Specific Constraints *****

Variables Obj;

Positive Variables
  purA           Purity of A
  recA           Recovery of A;

Scalars
  A_recoverySpec Recovery spec           /0.6/
  A_puritySpec   Purity spec             /0.551/;

Parameter
  ECost(u)       Equipment operating cost        /1 = 1.5, 2 = 1.0/;

Equations
  EqTotalFeed    Total feed flowrate of 1
  EqPurity_A     Purity calculations for A
  EqRecovery_A   Recovery calculations for A
  EqObj1         Minimize cost;

EqTotalFeed..
  SUM(FeedStr, F(FeedStr)) =e= 1;

EqPurity_A..
  purA*SUM(u, SUM(s$OutletV(u,s), F(s))) =e= SUM(u, SUM(s$OutletV(u,s), Fc(s,'A')));

EqRecovery_A..
  recA =e= SUM(u, SUM(s$OutletV(u,s), Fc(s,'A')))/FeedFlow('A');

EqObj1..
  Obj =e= SUM(u, SUM(s$Inlet(u,s), F(s)*ECost(u))) - 100*purA;

***** Initialization *****

* Assume a 50%/50% split of the feed between the two units

Fc.l(FeedStr,c) = FeedFlow(c)/2;
F.l(FeedStr) = SUM(c, Fc.l(FeedStr,c));
x.l(FeedStr,c) = Fc.l(FeedStr,c)/F.l(FeedStr);

* Intermediates for Newton step
Scalar g, dg, fd, xA;

* Vapor fraction
Scalar v, dv, iter;

* Feed composition
Parameter z(c);

* Initial vapor composition
Parameter vInit(u)       /1 = 0.5, 2 = 0.95/;

* Initialize each unit by performing VLE calculations
* Note: This initialization assumes two components, "A" and "B"
loop(u,
  v = vInit(u);
  iter = 0;
  g = 100;
  z(c) = SUM(s$Inlet(u,s), x.l(s,c));
  fd = SUM(s$Inlet(u,s), F.l(s));
  display z;
  while(abs(g) > 1E-8 AND iter < 20,
    xA = (1 - K(u,'B'))*z('B')/(( v*K(u,'B') + 1 - v)*( K(u,'A') - 1 ));
    g = v*K(u,'A')*xA + (1-v)*xA - z('A');
*   dg = K(u,'A') + (K(u,'B')-1)*xA*(1-v)/(1-v+K(u,'B')*v) - xA;
    dg = xA*(K(u,'A')- 1 + (1 - v + v*K(u,'A'))*(1 - K(u,'B'))/(1 - v + v*K(u,'B')) );
    dv = -g/dg;
    v = v + dv;
    v = max(min(v,1),0);
    iter = iter + 1;
    display g, dg, dv, v;
  );

  display xA, v, iter;

  loop(s$OutletL(u,s),
    F.l(s) = fd*(1-v);
    x.l(s,'A') = xA;
    x.l(s,'B') = 1-xA;
    Fc.l(s,c) = F.l(s)*x.l(s,c);
  );

  loop(s$OutletV(u,s),
    F.l(s) = fd*v;
    x.l(s,'A') = K(u,'A')*xA;
    x.l(s,'B') = 1 - x.l(s,'A');
    Fc.l(s,c) = F.l(s)*x.l(s,c);
  );

);


***** Bounds *****
F.lo(s) = 0;
x.lo(s,c) = 0.001;

purA.lo = A_puritySpec;
recA.lo = A_recoverySpec;

x.up(s,c) = 1;

x.fx(FeedStr,c) = FeedFlow(c)/SUM(c2, FeedFlow(c2));

option NLP = GAMSCHK;

Model Example1 /All/;

*option iterlim = 0;

Solve Example1 using NLP minimizing Obj;

option NLP = CONVERTD;
Example1.optfile = 1;

Solve Example1 using NLP minimizing Obj




