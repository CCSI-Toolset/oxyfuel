Pvap.l(Str, Comp) = exp(AntConst('1',Comp)+ AntConst('2',Comp)/T.l(Str)
                 + AntConst('5',Comp)*log(T.l(Str)) +  AntConst('6',Comp)*T.l(Str)**AntConst('7',Comp));

Pb.l(Str) = sum(Comp, Xc.l(Str,Comp)*Pvap.l(Str,Comp));
Pd.l(Str) = 1/sum(Comp, Xc.l(Str,Comp)/Pvap.l(Str,Comp));
