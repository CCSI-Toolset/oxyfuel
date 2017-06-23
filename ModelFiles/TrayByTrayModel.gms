**************************************************
* Tray-by-tray distillation modeling/optimization
**************************************************

Sets
  Trays                  Available distillation trays            /T1*T350/
  TrayLiqStr(Str)        Liquid streams for trays                /L1*L350/
  TrayVapStr(Str)        Vapor streams for trays                 /V1*V350/
  TrayStr(Str)           Streams for trays                       / /
  TrayCasc(Ed1,Trays)    Tray-cascade pairing                    / /
  ActvT(Trays)           Active distillation trays               / /
  InTrL(Trays,Str)       Inlet liquid streams for trays          / /
  OutTrL(Trays,Str)      Outlet liquid streams for trays         / /
  InTrV(Trays,Str)       Inlet vapor streams for trays           / /
  OutTrV(Trays,Str)      Outlet liquid streams for trays         / /
  BotTrays(Ed1,Trays)    Trays at bottoms of cascades            / /
  TopTrays(Ed1,Trays)    Trays at top of cascades                / /
  TopTraysOld(Ed1,Trays) Used for initialization adjustment      / /;

Positive Variables
  Kt(Trays,Comp)
  NumTrays(Ed1)
  eff(Trays);

Parameters
  NTrays                 Number of available trays               /350/
  CascSize(Ed1)          Cascade size (no. trays)                / /
  count                  counter - allocated trays               /0/;

*** Tray models ***

Equations
  EqTrayMoleBal(Trays,Comp)
  EqTrayEnrgBal(Trays)
  EqTrayKCEOS(Trays,Comp,Str,Str2)
  EqTrayXcEquilib(Trays,Comp,Str,Str2)
  EqTrayTEquilib(Trays,Str,Str2)
  EqTrayPEquilib(Trays,Str,Str2)
  EqTrayPDrop(Trays,Str,Str2)
  EqTraySumFrac(Trays,Str,Str2)
  EqTrayConstPres(Ed1,Str,Str2);

EqTrayMoleBal(Trays,Comp)..
  SUM(Str$(InTrL(Trays,Str) OR InTrV(Trays,Str)), Fc(Str,Comp)) =e=
    SUM(Str2$(OutTrL(Trays,Str2) OR OutTrV(Trays,Str2)), Fc(Str2,Comp));

EqTrayEnrgBal(Trays)..
  SUM(Str$(InTrL(Trays,Str) OR InTrV(Trays,Str)), F(Str)*H(Str)) =e=
    SUM(Str2$(OutTrL(Trays,Str2) OR OutTrV(Trays,Str2)), F(Str2)*H(Str2));

*ReturnHere
EqTrayKCEOS(Trays,Comp,Str,Str2)$(OutTrV(Trays,Str) AND OutTrL(Trays,Str2))..
  Kt(Trays,Comp)*phiEOS(Str,Comp) =e= phiEOS(Str2,Comp);

EqTrayXcEquilib(Trays,Comp,Str,Str2)$(OutTrV(Trays,Str) AND OutTrL(Trays,Str2))..
  Xc(Str,Comp) =e= Kt(Trays,Comp)*Xc(Str2,Comp);

EqTrayTEquilib(Trays,Str,Str2)$(OutTrV(Trays,Str) AND OutTrL(Trays,Str2))..
  Tscaled(Str) =e= Tscaled(Str2);

EqTrayPEquilib(Trays,Str,Str2)$(OutTrV(Trays,Str) AND OutTrL(Trays,Str2))..
  P(Str) =e= P(Str2);

EqTrayPDrop(Trays,Str,Str2)$(OutTrV(Trays,Str) AND InTrV(Trays,Str2) AND ActvT(Trays))..
  P(Str) =e= P(Str2) - eff(Trays)*SUM(Ed1$TrayCasc(Ed1,Trays), PresDropCnst(Ed1));

EqTrayConstPres(Ed1,Str,Str2)$(InVEd1(Ed1,Str) AND InLEd1(Ed1,Str2))..
  P(Str) =e= P(Str2) + NumTrays(Ed1)*PresDropCnst(Ed1);

EqTraySumFrac(Trays,Str,Str2)$(OutTrV(Trays,Str) AND OutTrL(Trays,Str2))..
  SUM(Comp, Xc(Str,Comp)) =e= SUM(Comp, Xc(Str2,Comp));

Model TrayEqns /EqTrayMoleBal, EqTrayEnrgBal, EqTrayKCEOS, EqTrayXcEquilib, EqTrayTEquilib, EqTrayPEquilib, EqTrayPDrop, EqTrayConstPres/;
* EqTraySumFrac
