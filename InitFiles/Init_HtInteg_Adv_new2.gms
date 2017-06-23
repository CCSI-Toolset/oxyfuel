loop(StrBase$PStr(StrBase),
  QAh.l(StrBase) = SUM(HtEx$CoolHtEx(HtEx), FCp.l(HtEx)*smmax(Tin.l(HtEx) - TPcnd.l(StrBase)) - FCp.l(HtEx)*smmax( Tout.l(HtEx) - TPcnd.l(StrBase) ) )
            +  SUM(TCond, FCp_cond.l(Tcond)*smmax(Tin_c.l(Tcond) - TPcnd.l(StrBase)) - FCp_cond.l(Tcond)*smmax( Tout_c.l(Tcond) - TPcnd.l(StrBase) ) );

  QAc.l(StrBase) = SUM(GnrlE$(HEqp.l(GnrlE) AND NOT InactiveGnrlE(GnrlE)), FCp.l(GnrlE)*smmax( Tout.l(GnrlE) - TPcnd.l(StrBase) + HRAT) - FCp.l(GnrlE)*smmax(Tin.l(GnrlE) - TPcnd.l(StrBase) + HRAT ) );

);

Qs.l = smax(PStr, QAc.l(PStr) - QAh.l(PStr));

Qw.l = Qs.l + SUM(HtEx, Qout.l(HtEx) - Qin.l(HtEx)) + SUM(TCond, QTCond.l(TCond)) - SUM(PReb, Qin.l(PReb));
