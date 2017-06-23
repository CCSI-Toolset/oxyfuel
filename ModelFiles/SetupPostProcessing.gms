************************************
*** Define Post Processing Macro ***
************************************

Variables liqPenAnlys(Str), vapPenAnlys(Str);

Variables Qin_tot, Qout_tot;

$macro ChckCmpl(x) liqPenAnlys.l(Str)$FlashLiq(Str)= F.l(Str)*sL.l(Str); \
         vapPenAnlys.l(Str)$FlashVap(Str) = F.l(Str)*sV.l(Str); \
         Qin_tot.l = SUM(GnrlE, Qin.l(GnrlE)); \
         Qout_tot.l = SUM(GnrlE, Qout.l(GnrlE)); \
         display liqPenAnlys.l, vapPenAnlys.l, Qin_tot.l, Qout_tot.l;
