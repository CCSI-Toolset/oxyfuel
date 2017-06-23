* Empty RelaxThermo
RelaxThermo(Str) = no;

loop(Str$(FlashVap(Str) AND RealStr(Str)),
  if(P.l(Str) > 1.01*Pd.l(Str) AND Pb.l(Str) - Pd.l(Str) > 0,
    RelaxThermo(Str) = yes;
  );
);

loop(Str$(FlashLiq(Str) AND RealStr(Str)),
  if(P.l(Str) < 0.99*Pb.l(Str) AND Pb.l(Str) - Pd.l(Str) > 0,
    RelaxThermo(Str) = yes;
  );
);

display RelaxThermo;
