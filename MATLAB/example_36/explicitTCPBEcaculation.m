function [xn1,xn2,xn3] = explicitTCPBEcaculation(beta,np,prFloc,prTime, ...
    xnc,fdia,xn1,xn2,xn3,shear,wsp,wsf)

brk_s=zeros(np,1);

% ------------------------------------------------------------------------
% Previous Scalar Concentration
ppxn1=xn1;
ppxn2=xn2;
ppxn3=xn3;

for k=1:np
    xnc(k)=xn3(k)/xn2(k);
    fdia(k)=xnc(k)^(1.0d0/prFloc.frac_df)*prFloc.pdia;

    beta.bm1(k)=2.0D0*prFloc.bolz*prFloc.temp/(3.0D0*prFloc.vmu)*4.0d0;
    beta.sh1(k)=1.0D0/6.0D0*shear(k)*(2.0D0*prFloc.pdia)^3.0D0;
    beta.pbe1(k)=prFloc.alpha_pbe1*(beta.bm1(k)+beta.sh1(k));

    beta.bm2(k)=2.0D0*prFloc.bolz*prFloc.temp/(3.0D0*prFloc.vmu)...
        *(1.0D0/prFloc.pdia+1.0D0/fdia(k))*(prFloc.pdia+fdia(k));
    beta.ds2(k)=3.1416D0/4.0D0*(fdia(k)+prFloc.pdia)^2.0D0...
        *abs(wsf(k)-wsp(k));
    beta.sh2(k)=1.0D0/6.0D0*shear(k)*(prFloc.pdia+fdia(k))^3.0D0;
    beta.pbe2(k)=prFloc.alpha_pbe2*(beta.bm2(k)+beta.ds2(k)+beta.sh2(k));

    beta.bm3(k)=2.0D0*prFloc.bolz*prFloc.temp/(3.0D0*prFloc.vmu)*4.0d0;
    beta.sh3(k)=1.0D0/6.0D0*shear(k)*(fdia(k)+fdia(k))^3.0D0;
    beta.pbe3(k)=prFloc.alpha_pbe3*(beta.bm3(k)+beta.ds3(k)+beta.sh3(k));

    brk_s(k)=prFloc.es*shear(k)*(fdia(k)/prFloc.pdia-1.0d0)^ ...
        (3.0D0-prFloc.frac_df) ...
        *(prFloc.vmu*shear(k)*fdia(k)^2.0D0*prFloc.ffy)^1.00d0;

    % return arguments
    xn1(k)=ppxn1(k)+prTime.dt*(-0.5D0*beta.pbe1(k)*ppxn1(k)^2*...
        (xnc(k)/(xnc(k)-1.0d0))-beta.pbe2(k)*ppxn1(k)*ppxn2(k)+...
        prFloc.brk_f*xnc(k)*brk_s(k)*ppxn2(k));
    xn2(k)=ppxn2(k)+prTime.dt*(0.5D0*beta.pbe1(k)*ppxn1(k)^2*...
        (1.0d0/(xnc(k)-1.0d0))-0.5D0*beta.pbe3(k)*ppxn2(k)^2+...
        brk_s(k)*ppxn2(k));
    xn3(k)=-xn1(k)+ppxn1(k)+ppxn3(k);
end

end

