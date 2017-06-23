*************************
***** InitThrmE.gms *****
*************************

* Assume K has been calculated in the thermo initialization script

display beta.l, CalcThrmE;

loop(ThrmE$CalcThrmE(ThrmE),
  beta.l(ThrmE) = SUM(Comp, SUM(Str$OutVThrmE(ThrmE,Str), Xc.l(Str,Comp))/
         max( K.l(ThrmE,Comp)*SUM(Str$OutLThrmE(ThrmE,Str), Xc.l(Str,Comp)),1E-3))/CARD(Comp);
);

display beta.l;
