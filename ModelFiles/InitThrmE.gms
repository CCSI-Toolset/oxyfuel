*************************
***** InitThrmE.gms *****
*************************

* This script initializes each
* ThrmE unit by calculating the VLE output by iteration

scalar countT;
Parameter oldY(Comp), oldX(Comp);

loop(ThrmE,

   countT = 0;
   while(countT < 5,
     oldY(Comp) = Xc.l(Str,Comp)$OutVThrmE(ThrmE,Str);
     oldX(Comp) = Xc.l(Str,Comp)$OutLThrmE(ThrmE,Str);

     if(ThermoSwitch eq 1,
       Pvap.l(Str, Comp)$(OutLThrmE(ThrmE,Str) OR OutVThrmE(ThrmE,Str)) =
         exp(AntConst('1',Comp)+ AntConst('2',Comp)/T.l(Str)
         + AntConst('5',Comp)*log(T.l(Str)) +  AntConst('6',Comp)*T.l(Str)**AntConst('7',Comp));

       K.l(ThrmE,Comp) = SUM(Str$OutVThrmE(ThrmE,Str), Pvap.l(Str,Comp)/P.l(Str));

     else
       K.l(ThrmE,Comp) = 1;
     );

     Xc.l(Str,Comp)$OutVThrmE(ThrmE,Str) = K.l(ThrmE,Comp)*Xc.l(Str2,Comp)$OutLThrmE(ThrmE,Str2);


   );
);