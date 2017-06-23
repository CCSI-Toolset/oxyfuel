FcTarget(Str,Comp)$InitStr1(Str) = Fc.l(Str,Comp);
Ttarget(Str)$(InitStr1(Str) OR InitStr2(Str)) = Tscaled.l(Str);
FTarget(Str)$InitStr2(Str) = F.l(Str);
PTarget(Str)$(InitStr1(Str) OR InitStr2(Str)) = P.l(Str);
XcTarget(Str,Comp)$InitStr1(Str) = Xc.l(Str,Comp);

FPen.l = SUM(Str$InitStr2(Str), (F.l(Str) - FTarget(Str))*(F.l(Str) - FTarget(Str)));
FcPen.l = SUM(Comp, SUM(Str$InitStr1(Str), (Fc.l(Str,Comp) - FcTarget(Str,Comp))*(Fc.l(Str,Comp) - FcTarget(Str,Comp))));
TPen.l = SUM(Str$(InitStr1(Str) OR InitStr2(Str)), (Ttarget(Str) - Tscaled.l(Str))*(Ttarget(Str) - Tscaled.l(Str)) );
PPen.l = SUM(Str$(InitStr1(Str) OR InitStr2(Str)), (PTarget(Str) - P.l(Str))*(PTarget(Str) - P.l(Str)));
XcPen.l = SUM(Comp, SUM(Str$InitStr1(Str), (XcTarget(Str,Comp) - Xc.l(Str,Comp))*(XcTarget(Str,Comp) - Xc.l(Str,Comp)) ));