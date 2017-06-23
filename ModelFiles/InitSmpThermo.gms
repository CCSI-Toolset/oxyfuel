*****************************
***** InitSmpThremo.gms *****
*****************************

* This file initializes the simple thermodynamic variables

Pvap.l(Str, Comp) = exp(AntConst('1',Comp)+ AntConst('2',Comp)/T.l(Str)
                 + AntConst('5',Comp)*log(T.l(Str)) +  AntConst('6',Comp)*T.l(Str)**AntConst('7',Comp));

Pb.l(Str) = sum(Comp, Xc.l(Str,Comp)*Pvap.l(Str,Comp));
Pd.l(Str) = 1/sum(Comp, Xc.l(Str,Comp)/Pvap.l(Str,Comp));

K.l(ThrmE,Comp) = SUM(Str$OutVThrmE(ThrmE,Str), Pvap.l(Str,Comp)/P.l(Str));

beta.l(ThrmE) = 1;

$include ./IncludeFiles/Init_HIG.gms

sV.l(FlashVap) = 1;
sL.l(FlashLiq) = 1;