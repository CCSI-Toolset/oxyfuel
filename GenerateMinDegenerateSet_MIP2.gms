************************************************************
***** MILP for Irreducible Set of Degenerate Equations *****
************************************************************

* Created by Alex Dowling at Carnegie Mellon University
* awdowlin@andrew.cmu.edu
* Last Modified: June 10th, 2014

$onempty

$set matout "'degenSol.gdx', x, y, Z, MIPstat ";

Sets
  i              Equations                       /e1*e100000/
  j              Variables                       /x1*x100000/
  actE(i)        Equations in active Jacobian
  SE(i)          Suspect equation
;

* Note: Increasing the number of equations and variables declared in GAMS (above) increases run time.
* Currently 100,000 variables and equations are supported. To increase this number, simply increase the
* size of i & j.

Parameter
  A(i,j) Jacobian
  MIPstat
;

$gdxin jacobian.gdx
$load A
$gdxin

$gdxin degenData.gdx
$load actE
$load SE
$gdxin

display actE, SE;

*stop

Scalar
  bigM   Big-M constant  /10/
;

Variable
  x(i)   Eigenvector
  z      Objective function value;

Binary Variable
  y(i)   Vector element toggle;

Equations
  EqDegenecary(j)
  xlow(i)
  xup(i)
  EqObj ;

EqDegenecary(j)..
  SUM(i$actE(i), A(i,j)*x(i)) =e= 0;

xlow(i)$actE(i)..
  -bigM*y(i) =l= x(i);

xup(i)$actE(i)..
  x(i) =l= bigM*y(i);

EqObj..
  z =e= SUM(i$actE(i), y(i));

Model GenerateMinDegenerateSet /All/;

*option MINLP = SBB;
*option NLP = SNOPT;
*option MINLP = DICOPT;

*option MIP = GUROBI;
option MIP = CPLEX;

GenerateMinDegenerateSet.optca = 0.9;
GenerateMinDegenerateSet.optcr = 0.01;

x.l(i) = 1;
y.l(i) = 1;

y.fx(SE) = 1;
x.fx(SE) = 1;

*GenerateMinDegenerateSet.optfile = 1;
Solve GenerateMinDegenerateSet minimizing Z using MIP;

MIPstat = GenerateMinDegenerateSet.modelStat;

execute_unload %matout%;

