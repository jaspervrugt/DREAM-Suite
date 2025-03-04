function [xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,wsp,wsf,ers_np,dep_nt]...
    = implicitSedTransCalculation...
    (ae2,ae3,np,tau_b,delz,prFloc,prTime,prBot,...
    xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,tde1,wsp,wsf)

% ------------------------------------------------------------------------
% Initialization for sediment transport
change_max=1.0d0;
j=0;

pxn1=xn1;
pxn2=xn2;
pxn3=xn3;

% ------------------------------------------------------------------------
% deposition/erosion (sink/source) calculation - Le Hir (2011, csr)
if (tau_b>prBot.tau_c)
    ers_np=prBot.ers_m/(1.0d0/6.0d0*3.1416d0*(prFloc.pdia^3)*prFloc.par_den)...
        *(tau_b/prBot.tau_c-1.0d0);
else
    ers_np=0.0d0;
end
dep_np=wsp(1)*xn1(1);
dep_nf=wsf(1)*xn2(1);
dep_nt=wsf(1)*xn3(1);
ers_np=ers_np+dep_np;

while change_max > 1.0d-10
    j=j+1;
    cxn3=xn3; % Just for checking numerical solution
    % ------------------------------------------------------------------------
    % update mass and volume concentration
    % update settlinprFloc.g velocity: winterwerp (2004) and r-z eqn (burprFloc.ger, 2008)
    for k=1:np
    	xnc(k)=xn3(k)/xn2(k);
     	fdia(k)=xnc(k)^(1.0d0/prFloc.frac_df)*prFloc.pdia;
    	xnt(k)=xn1(k)+xn3(k);
    	conc(k)=xnt(k)*prFloc.par_den*prFloc.pri_vol;
    	floc_den(k)=prFloc.wat_den+(prFloc.par_den-prFloc.wat_den)...
         	*(prFloc.pdia/fdia(k))^(3.0d0-prFloc.frac_df);
    	wsf(k)=-0.055556d0*prFloc.g*(floc_den(k)/prFloc.wat_den-1.0d0)...
         	*fdia(k)^2/(prFloc.vmu*1.0d-3);
    	wsp(k)=-0.055556d0*prFloc.g*(prFloc.par_den/prFloc.wat_den-1.0d0)*prFloc.pdia^2/(prFloc.vmu*1.0d-3);
    end

    % ------------------------------------------------------------------------
    % (note) diffusion coeff. for 2nd order approximation wrt time
    for k=1:np
        if (k==1)
    	    ae2(k)=prTime.dt*(tde1(k+1)+tde1(k))/(2.0d0*delz(k)^2);
  	  elseif (k==np)
  	    ae3(k)=prTime.dt*(tde1(k)+tde1(k-1))/(2.0d0*delz(k-1)^2);
        else
    	    ae2(k)=prTime.dt*(tde1(k+1)+tde1(k))/(2.0d0*delz(k)^2);
    	    ae3(k)=prTime.dt*(tde1(k)+tde1(k-1))/(2.0d0*delz(k-1)^2);
        end
    end


    % scalar (sediment concentration) update inside nodal system
    % ------------------------------------------------------------------------
    % node calculation at the bottom boundary (open)
    for k=1:np
        if (k==1)
      	  xn1(k)=(pxn1(k)-...
           	prTime.dt/delz(k)*(wsp(k+1)*xn1(k+1)-ers_np)+ae2(k)*xn1(k+1))/...
         	(1.0d0+ae2(k));
    	  xn2(k)=(pxn2(k)-...
           	prTime.dt/delz(k)*(wsf(k+1)*xn2(k+1)-dep_nf)+ae2(k)*xn2(k+1))/...
         	(1.0d0+ae2(k));
    	  xn3(k)=(pxn3(k)-...
           	prTime.dt/delz(k)*(wsf(k+1)*xn3(k+1)-dep_nt)+ae2(k)*xn3(k+1))/...
         	(1.0d0+ae2(k));
          % ------------------------------------------------------------------------
          % node calculation at the top boundary (closed)
        elseif (k==np)
      	  xn1(k)=(pxn1(k)+ae3(k)*xn1(k-1))/...
           	(1.0d0-prTime.dt/delz(k)*wsp(k)+ae3(k));
    	  xn2(k)=(pxn2(k)+ae3(k)*xn2(k-1))/...
           	(1.0d0-prTime.dt/delz(k)*wsf(k)+ae3(k));
    	  xn3(k)=(pxn3(k)+ae3(k)*xn3(k-1))/...
           	(1.0d0-prTime.dt/delz(k)*wsf(k)+ae3(k));
          % ------------------------------------------------------------------------
          % node calculation inside boundary
        else
      	  xn1(k)=(pxn1(k)-prTime.dt/delz(k)*wsp(k+1)*xn1(k+1)...
           	+ae2(k)*xn1(k+1)+ae3(k)*xn1(k-1))/...
         	(1.0d0-prTime.dt/delz(k)*wsp(k)+ae2(k)+ae3(k));
    	  xn2(k)=(pxn2(k)-prTime.dt/delz(k)*wsf(k+1)*xn2(k+1)...
           	+ae2(k)*xn2(k+1)+ae3(k)*xn2(k-1))/...
         	(1.0d0-prTime.dt/delz(k)*wsf(k)+ae2(k)+ae3(k));
    	  xn3(k)=(pxn3(k)-prTime.dt/delz(k)*wsf(k+1)*xn3(k+1)...
           	+ae2(k)*xn3(k+1)+ae3(k)*xn3(k-1))/...
         	(1.0d0-prTime.dt/delz(k)*wsf(k)+ae2(k)+ae3(k));
        end
    end

    % ------------------------------------------------------------------------
    % check converprFloc.gence
    change=abs((cxn3-xn3)./xn3);
    change_max=max(change);

    % ------------------------------------------------------------------------
    % end of prFloc.gauss-siedel iteration
end

end
