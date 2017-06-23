*******************************
***** FinalProcessing.gms *****
*******************************

* This file includes a few post processing calculations:
*  - Calculate combined stream (ThrmE outlet vapor & liquid) properties
*      * Flow
*      * Composition
*      * Vapor fraction
*  - Vapor pressure calculations - useful for phase checks/comparisons
*  - Entropy balance checks for ThrmE and cascades

Variables
  F_ss(ThrmE)
  T_ss(ThrmE)
  P_ss(ThrmE)
  Fc_ss(ThrmE,Comp)
  Xc_ss(ThrmE,Comp)
  VapFrac(ThrmE);

F_ss.l(ThrmE) = SUM(Str$(OutVThrmE(ThrmE,Str) OR OutLThrmE(ThrmE,Str)), F.l(Str));
T_ss.l(ThrmE) = SUM(Str$OutVThrmE(ThrmE,Str), T.l(Str));
P_ss.l(ThrmE) = SUM(Str$OutVThrmE(ThrmE,Str), P.l(Str));
Fc_ss.l(ThrmE,Comp) = SUM(Str$(OutVThrmE(ThrmE,Str) OR OutLThrmE(ThrmE,Str)), Fc.l(Str,Comp));
Xc_ss.l(ThrmE,Comp) = Fc_ss.l(ThrmE,Comp)/max(F_ss.l(ThrmE),1E-10);
VapFrac.l(ThrmE) = SUM(Str$OutVThrmE(ThrmE,Str), F.l(Str))/max(F_ss.l(ThrmE),1E-10);

display F_ss.l, T_ss.l, P_ss.l, Fc_ss.l, Xc_ss.l, VapFrac.l;

Pvap.l(Str, Comp) = exp(AntConst('1',Comp)+ AntConst('2',Comp)/T.l(Str)
+ AntConst('5',Comp)*log(T.l(Str)) +  AntConst('6',Comp)*T.l(Str)**AntConst('7',Comp));

Pb.l(Str) = sum(Comp, Xc.l(Str,Comp)*Pvap.l(Str,Comp));
Pd.l(Str) = 1/sum(Comp, Xc.l(Str,Comp)/Pvap.l(Str,Comp));

display Pvap.l, Pb.l, Pd.l;

LogS1.l(fEOSStr) = max((ZEOS.l(fEOSStr) - bbEOS.l(fEOSStr))/(ZEOS.l(fEOSStr)),1E-10);
LogS2.l(fEOSStr) = max(ZEOS.l(fEOSStr)*Pref/P.l(fEOSStr),1E-10);

S.l(fEOSStr) = SIG.l(fEOSStr) + Rsi*log(Pref/P.l(fEOSStr)) + 100*(R*log(LogS1.l(fEOSStr)) +
         R*log(LogS2.l(fEOSStr)) + TdadT.l(fEOSStr)/(bmEOS.l(fEOSStr)*T.l(fEOSStr))
         *log(LogS3.l(fEOSStr))) ;

Variables
  Sin(ThrmE),
  Sout(ThrmE),
  Sgen(ThrmE),
  Sin_c(Ed1),
  Sout_c(Ed1),
  Sgen_c(Ed1);

loop(ThrmE,
  Sin.l(ThrmE) = SUM(Str$InThrmE(ThrmE,Str), F.l(Str)*S.l(Str));
  Sout.l(ThrmE) = SUM(Str$(OutVThrmE(ThrmE, Str) OR OutLThrmE(ThrmE, Str)), F.l(Str)*S.l(Str));
  Sgen.l(ThrmE) = Sout.l(ThrmE) - Sin.l(ThrmE);
);

loop(Ed1,
  Sin_c.l(Ed1) = SUM(Str$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str)), F.l(Str)*S.l(Str));
  Sout_c.l(Ed1) = SUM(Str$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)), F.l(Str)*S.l(Str));
  Sgen_c.l(Ed1) = Sout_c.l(Ed1) - Sin_c.l(Ed1);
);

display LogS1.l, LogS2.l, SIG.l, S.l, Sin.l, Sout.l, Sgen.l, Sin_c.l, Sout_c.l, Sgen_c.l;