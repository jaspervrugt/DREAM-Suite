%% Connector to DREAM-MCMC
function [ ustar,tau_b,tde1,te,shear ] =  ...
    updataFlowfield (time,prFloc,prGrid,prFlow,tde1,te,shear)

np = prGrid.np;

% ------------------------------------------------------------------------
% bottom shear (ustar) and shear rate profile
	ustar=abs(0.036d0+0.000d0 ...
        *sin(2.0d0*3.1415d0/12.4d0*(time/3.6d3)+0.5d0) ...
        -0.020d0 ...
        *cos(2.0d0*(2.0d0*3.1415d0/12.4d0*(time/3.600d3)+0.5d0)));
	tau_b=prFloc.wat_den*ustar^2;
% ------------------------------------------------------------------------
% update eddy viscosity and turbulence-related parameter

if (np==1) % Depth-averaged Model
      k = 1;
	  tde1(k)=prFlow.v_karm*ustar*double(k)*(prGrid.dz/2.0)*(1.0d0-double(k)*(prGrid.dz/2)/prGrid.t_height);
	  te(k)=ustar^3/(prFlow.v_karm*(prGrid.dz/2)*double(k));
	  shear(k)=sqrt(te(k)/(prFloc.vmu*1.0d-3));
else
    for k=1:np
	  tde1(k)=prFlow.v_karm*ustar*double(k)*prGrid.dz*(1.0d0-double(k)*prGrid.dz/prGrid.t_height);
	  te(k)=ustar^3/(prFlow.v_karm*prGrid.dz*double(k));
	  shear(k)=sqrt(te(k)/(prFloc.vmu*1.0d-3));
    end
end    


