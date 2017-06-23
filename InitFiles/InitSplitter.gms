loop(Sptr,
  P.l(Str)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), P.l(Str2));
  Tscaled.l(Str)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), Tscaled.l(Str2));
  H.l(Str)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), H.l(Str2));
  phiEOS.l(Str,Comp)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), phiEOS.l(Str2,Comp));
  ZEOS.l(Str)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), ZEOS.l(Str2));
  Pvap.l(Str,Comp)$OutSptr(Sptr,Str) = SUM(Str2$InSptr(Sptr,Str2), Pvap.l(Str2,Comp));
);