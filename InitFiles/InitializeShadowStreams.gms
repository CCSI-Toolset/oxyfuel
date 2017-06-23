***************************************
***** InitializeShadowStreams.gms *****
***************************************

* Use Newton's method and Antiones eqn for dew/bubble point calculations

* Difficult part.. shadow streams for phase stability analysis
* This require Newton's method to calculate Tdew or Tbub

PrmShd(Str) = no;

loop(Str$PhaseStability(Str),
  loop(VapShd,
    loop(LiqShd$ComMap(Str,VapShd,LiqShd),
      if(VapStr(Str) OR FlashVap(Str),
        PrmShd(VapShd) = yes;
        Xc.l(VapShd,comp) = Xc.l(Str,Comp);
      elseif LiqStr(Str) OR FlashLiq(Str),
        PrmShd(LiqShd) = yes;
        Xc.l(LiqShd,Comp) = Xc.l(Str,Comp);
      );
      Tscaled.l(VapShd) = Tscaled.l(Str);
      Tscaled.l(LiqShd) = Tscaled.l(Str);
      P.l(VapShd) = P.l(Str);
      P.l(LiqShd) = P.l(Str);
    );
  );
);

* Solve for temperature in shadow streams involved in PhaseStability
for(iter = 0 to 10,
  Pvap.l(PrmShd,Comp) = exp(AntConst('A',Comp) - AntConst('B',Comp)/(Tscaled.l(PrmShd)*Tref + AntConst('C',Comp)))/750.061683;
  dPvapdT(PrmShd,Comp)= AntConst('B',Comp)*Pvap.l(PrmShd,Comp)/Power(Tscaled.l(PrmShd)*Tref + AntConst('C',Comp), 2);

    Pdtemp(PrmShd) = 1/SUM(Comp, Xc.l(PrmShd,Comp)/Pvap.l(PrmShd,Comp));
    Pbtemp(PrmShd) = SUM(Comp, Xc.l(PrmShd,Comp)*Pvap.l(PrmShd,Comp));

    loop(PrmShd,
      if( VapShd(PrmShd),
        dTs(PrmShd) = -(Pdtemp(PrmShd) - P.l(PrmShd))
          /(SUM(Comp, Xc.l(PrmShd,Comp)/Power(Pvap.l(PrmShd,Comp),2)*dPvapdT(PrmShd,Comp))*Power(Pdtemp(PrmShd),2));
      elseif LiqShd(PrmShd),
        dTs(PrmShd) = - (Pbtemp(PrmShd) - P.l(PrmShd))/SUM(Comp, Xc.l(PrmShd,Comp)*dPvapdT(PrmShd,Comp));
      );
    );
    Tscaled.l(PrmShd) = Tscaled.l(PrmShd) + dTs(PrmShd)/Tref;

    if(DebugMode > 2,
      display dTs;
    );
  );

loop(Str$PhaseStability(Str),
  loop(PrmShd,
    loop(Str2$(ComMap(Str,PrmShd,Str2) OR ComMap(Str,Str2,PrmShd)),
      Tscaled.l(Str2) = Tscaled.l(PrmShd);
    );
  );
);

* Easy part... calculate equilibrium compositions
* Assume DewPoint( ) and BubPoint( ) streams are already at their bubble or dew
* point temperature
loop(Str$(DewPoint(Str) OR BubPoint(Str)),
  loop(VapShd,
    loop(LiqShd$ComMap(Str,VapShd,LiqShd),
      Xc.l(VapShd,Comp) = Pvap.l(Str,Comp)*Xc.l(Str,Comp)/P.l(Str);
      Xc.l(LiqShd,Comp) = Xc.l(Str,Comp)*P.l(Str)/Pvap.l(Str,Comp);
      Tscaled.l(VapShd) = Tscaled.l(Str);
      Tscaled.l(LiqShd) = Tscaled.l(Str);
      P.l(VapShd) = P.l(Str);
      P.l(LiqShd) = P.l(Str);
    );
  );
);

* Another easy part... calculate the equilibrium compositions for PhaseStability
* shadow streams.
loop(Str$PhaseStability(Str),
  loop(PrmShd,

* Dew point calculation... initialize liquid shadow stream
    loop(Str2$ComMap(Str,PrmShd,Str2),
      Xc.l(Str2,Comp) = Xc.l(PrmShd,Comp)*P.l(PrmShd)/Pvap.l(PrmShd,Comp);
    );

* Bubble point calculation... initialize vapor shadow stream
    loop(Str2$ComMap(Str,Str2,PrmShd),
      Xc.l(Str2,Comp) = Pvap.l(PrmShd,Comp)*Xc.l(PrmShd,Comp)/P.l(PrmShd);
    );
  );
);
