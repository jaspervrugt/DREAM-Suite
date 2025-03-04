function [xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,wsf,wsp,pri_tmass, ...
    tde1,te,shear] = InitialConditions(time,prIni,prFloc,prGrid, ...
    prFlow,pri_tmass)
% Function defines the initial condition

np = prGrid.np;

% ------------------------------------------------------------------------
% Preallocate all the arrays
[xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den, ...
    wsp,wsf,tde1,te,shear] = deal(zeros(np,1));

ustar0 = abs(0.036d0+0.000d0 ...
    * sin(2.0d0*3.1415d0/12.4d0*(time/3.6d3)+0.5d0)...
    -0.020d0 ...
    * cos(2.0d0*(2.0d0*3.1415d0/12.4d0*(time/3.600d3)+0.5d0)));

for k=1:np
    % sediment and floc profilprFloc.es
    xn1(k)=prIni.pnconc*0.90d0;
    xn3(k)=prIni.pnconc*0.10d0;
    fdia(k)=prFloc.ffdia;
    xnc(k)=(fdia(k)/prFloc.pdia)^prFloc.frac_df;
    xn2(k)=xn3(k)/xnc(k);
    xnt(k)=xn1(k)+xn3(k);
    conc(k)=xnt(k)*prFloc.par_den*prFloc.pri_vol;
    floc_den(k)=prFloc.wat_den+(prFloc.par_den-prFloc.wat_den)...
     	*(prFloc.pdia/fdia(k))^(3.0d0-prFloc.frac_df);
    re_p=fdia(k)*-wsf(k)*prFloc.wat_den/prFloc.vmu;
    wsf(k)=-(prFloc.par_den-prFloc.wat_den)*prFloc.g * ...
        prFloc.pdia^(3-prFloc.frac_df)...
        *fdia(k)^(prFloc.frac_df-1)/(18.0d0*prFloc.vmu)...
        /(1+0.15*re_p^0.687);
    wsp(k)=-(prFloc.pdia^2)*prFloc.g*(prFloc.par_den - ...
        prFloc.wat_den)/(18.0d0*prFloc.vmu);
    pri_tmass=pri_tmass+xnt(k)*prGrid.dz;
    %  f_hind(k) = 1.0d0;
    % flow and turbulence profile
    tde1(k)=prFlow.v_karm*ustar0*double(k)*prGrid.dz*(1.0d0 - ...
        double(k)*prGrid.dz/prGrid.t_height);
    te(k)=ustar0^3/(prFlow.v_karm*prGrid.dz*double(k));
    shear(k)=sqrt(te(k)/(prFloc.vmu*1.0d-3));
end

end
