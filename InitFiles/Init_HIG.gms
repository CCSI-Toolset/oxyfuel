Tscaled.l(Str) = T.l(Str)/Tref;

HIG.l(Str) = 1E-3*SUM(J, Xc.l(Str,j)*(Tref**5 *CpIG('5',J)/5*(Tscaled.l(Str)**5 - 1) + Tref**4 *CpIG('4',J)/4*(Tscaled.l(Str)**4 - 1)
         + Tref**3 *CpIG('3',J)/3*(Tscaled.l(Str)**3 - 1) + Tref**2 *CpIG('2',J)/2*(Tscaled.l(Str)**2 - 1) + Tref*CpIG('1',J)*(Tscaled.l(Str) - 1)));
