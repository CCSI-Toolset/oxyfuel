*****************************
***** InitSmpThremo.gms *****
*****************************

* This file initializes the simple thermodynamic variables

*Tscaled.l(Str) = T.l(Str)/Tref;
T.l(Str) = Tscaled.l(Str)*Tref;

T.l(Str) = max(T.l(Str),1);

display T.l;

Pvap.l(Str, Comp) = exp(AntConst('A',Comp) - AntConst('B',Comp)/(Tscaled.l(Str)*Tref + AntConst('C',Comp)))/750.061683;

Pb.l(Str) = sum(Comp, Xc.l(Str,Comp)*Pvap.l(Str,Comp));
Pd.l(Str) = 1/sum(Comp, Xc.l(Str,Comp)/Pvap.l(Str,Comp));

K.l(ThrmE,Comp) = SUM(Str$OutVGnrlE(ThrmE,Str), Pvap.l(Str,Comp)/P.l(Str));

beta.l(ThrmE) = 1;

$include ../InitFiles/Init_HIG.gms

sV.l(FlashVap) = 1;
sL.l(FlashLiq) = 1;

***** TODO - Check this

loop(Str,
  if(P.l(Str) <= Pd.l(Str) + epsi,
    sL.l(Str) = 1;
  );
  if(P.l(Str) >= Pb.l(Str) - epsi,
    sV.l(Str) = 1;
  );
);

* Initialize enthalpy

H.l(Str)$(VapStr(Str) + FlashVap(Str)) = HIG.l(Str);

H.l(Str)$(LiqStr(Str) + FlashLiq(Str)) = HIG.l(Str) - SUM(Comp, Xc.l(Str,Comp)*1E-3*Rsi*(AntConst('B',Comp)/Power(Tscaled.l(Str)*Tref + AntConst('C',Comp),2) ));


H.l(Str)$(VapStr(Str) + FlashVap(Str)) = SUM(Comp, Xc.l(Str,Comp)*
         (HVapSurf('1',Comp) + HVapSurf('2',Comp)*P.l(Str) + HVapSurf('3',Comp)*Power(P.l(Str),2)
         + HVapSurf('4',Comp)*Tscaled.l(Str) + HVapSurf('5',Comp)*Power(Tscaled.l(Str),2)
         + HVapSurf('6',Comp)*Power(Tscaled.l(Str),3) + HVapSurf('7',Comp)*P.l(Str)*Tscaled.l(Str)
         + HVapSurf('8',Comp)*P.l(Str)*Power(Tscaled.l(Str),2)));

H.l(Str)$(LiqStr(Str) + FlashLiq(Str)) = SUM(Comp, Xc.l(Str,Comp)*(HLiqSurf('1',Comp)
         + HLiqSurf('4',Comp)*Tscaled.l(Str) + HLiqSurf('5',Comp)*Power(Tscaled.l(Str),2)
         + HLiqSurf('6',Comp)*Power(Tscaled.l(Str),3) + HLiqSurf('7',Comp)*P.l(Str)*Tscaled.l(Str)
         + HLiqSurf('8',Comp)*P.l(Str)*Power(Tscaled.l(Str),2)));
