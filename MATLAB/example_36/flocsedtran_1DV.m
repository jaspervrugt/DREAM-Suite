function floc_sim = flocsedtran_1DV(par,plugin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% 2-parameter flocullation model                                        %%
%%                                                                       %%
%%  SYNOPSIS: floc_sim = flocsedtran_1DV(par,plugin)                     %%
%% where                                                                 %%
%%  x           [input] 1xd vector of parameters                         %%
%%        M        Empirical erosion parameter          (kg/m2/s)        %%
%%        τ_c      Critical shear stress for erosion    (Pa)             %%
%%  plugin      [input] structure of model inputs                        %%
%%  simdata     [outpt] simulated temperature data                       %%
%%                                                                       %%
%%  LITERATURE                                                           %%
%%                                                                       %%
%%  © Written by B.J. Lee & modified by Jasper A. Vrugt, Aug. 2024       %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Unpack parameter values
prFloc = plugin.prFloc;
prFloc.alpha_pbe1 = par(1);
prFloc.alpha_pbe2 = par(1);
prFloc.alpha_pbe3 = par(1);
prFloc.frac_df = par(2);
np = plugin.prGrid.np;
% Initial states
ts_onemab = plugin.ts_onemab;
% initial condition (intial seeding)
sig_b= plugin.sig_b; time = plugin.time; 
% Execute Initial Conditions
[xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,wsf,wsp,pri_tmass,tde1, ...
    te,shear] = InitialConditions(time,plugin.prIni,prFloc, ...
    plugin.prGrid,plugin.prFlow,0);

for itm = 1:32400     % Time loop
    time = double(time+plugin.prTime.dt);
    % set the previous time value of bottom layer
    psig_b = sig_b;
    % Update Flow Fields
    [ustar,tau_b,tde1,te,shear] = ...
        updataFlowfield(time,prFloc,plugin.prGrid,plugin.prFlow, ...
        tde1,te,shear);
    % Explicit solver for tcpbe
    [xn1,xn2,xn3] = explicitTCPBEcaculation...
        (plugin.beta,np,prFloc,plugin.prTime,xnc,fdia,xn1,xn2,xn3, ...
        shear,wsp,wsf);
    % GS iteration for solving implicit 1st order sed transport equation
    if (np==1)
        [xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,wsp,wsf,ers_np,dep_nt]...
            = depthaveSedtran...
            (np,tau_b,plugin.delz,prFloc, ...
            plugin.prTime,plugin.prBot,...
            xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,tde1,wsp,wsf) ;
    else
        [xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,wsp,wsf,ers_np,dep_nt] ...
            = implicitSedTransCalculation...
            (plugin.ae2,plugin.ae3,np,tau_b,plugin.delz,prFloc, ...
            plugin.prTime,plugin.prBot,...
            xn1,xn2,xn3,xnt,conc,fdia,xnc,floc_den,tde1,wsp,wsf);
    end
    % Update thickness of bottom layer (bl)
    sig_b = psig_b + plugin.prTime.dt * plugin.prBot.fntc / ...
        prFloc.cgel*(-ers_np-dep_nt);
    % Post processing at each time
    if any(plugin.ts_obs >= time - plugin.prTime.dt & plugin.ts_obs < time)
        onemab_pfrac = xn1(1)/(xn1(1)+xn3(1));
        onemab_fdia = fdia(1)*1.0d3; % millimeter
        onemab_conc = conc(1);
        timeHour = (time)/3600;
        onemab_combined = cat...
            (2,timeHour,ustar,onemab_pfrac,onemab_fdia,onemab_conc);
        ts_onemab = cat(1,ts_onemab,onemab_combined);
    end

end % End of the time loop
ts_onemab
% Post processing at the end of time loop
ts_onemab = ts_onemab(2:end,3:end); % Removal the first row
floc_sim = ts_onemab(:);            % Return as a single vector

end
