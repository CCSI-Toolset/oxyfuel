*************************************
***** Degeneracy Hunter Example *****
*************************************

* Overspecified CSTR
* Created by Alex Dowling at Carnegie Mellon University
* awdowlin@andrew.cmu.edu
* Last Modified: June 10th, 2014

Scalars
  k              Rate constant                   /1.5/
  v              Feed and outlet velocity        /2.0/
  ;

Sets
  s              Streams                         /in, out/
  c              Components                      /A, B, C/
  ;

Variables
  F(s)           Total molar flowrate
  Fc(s,c)        Component molar flowrate
  xc(s,c)        Mole fraction
  Vol            CSTR volume
  X              CSTR conversion
  r              Reaction rate
  Ca             Concentration of A exiting reactor
  ;

Equations
  TotalStreamFlow(s)
  StreamComp(s,c)
  SumMoleFrac(s)
  CSTRDesignEqn
  ConvDef
  FormB
  RateLaw
  StreamConc
  Inerts
  ;

TotalStreamFlow(s)..
  F(s) =e= SUM(c, Fc(s,c));

StreamComp(s,c)..
  Fc(s,c) =e= F(s) * xc(s,c);

SumMoleFrac(s)..
  1 =e= sum(c, xc(s,c));

CSTRDesignEqn..
  -r * Vol =e= Fc('in','A') * X;

ConvDef..
  Fc('in','A') * X =e= Fc('in','A') - Fc('out','A');

FormB..
  Fc('out','B') - Fc('in','B') =e= Fc('in','A') - Fc('out','A');

RateLaw..
  r =e= -k * Ca;

StreamConc..
  Fc('out','A') =e= v * Ca;

Inerts..
  Fc('in','C') =e= Fc('out','C');

***** Bounds *****

Fc.lo(s,c) = 0;
xc.lo(s,c) = 0;
xc.up(s,c) = 1;
Vol.lo = 0;
Vol.up = 10;

***** Fix Variables *****
F.fx('in') = 1;
xc.fx('in','A') = 0.8;
xc.fx('in','C') = 0.2;


Model CSTRExample /all/

option NLP = CONOPT;
solve CSTRExample maximizing X using NLP;

option NLP = CONVERTD;
CSTRExample.optfile = 1;
solve CSTRExample maximizing X using NLP;

