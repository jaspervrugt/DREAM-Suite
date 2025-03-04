function plugin = model_properties
% Floc/Particle propertiflocParam.es

% Aggregation/Break-up Kinetic Constants
prFloc.alpha_pbe1=0.10D0;
prFloc.alpha_pbe2=0.10D0;
prFloc.alpha_pbe3=0.10D0;

% Correction Coefficient for Diff Settling-Mediated Flocculation
prFloc.es=2.0D-5;
prFloc.ffy=1.0D+10;
prFloc.brk_f=1.0D0;

% Constant for Fractal Theory - Check FlprFloc.esch et al (AICHE, 1999)
prFloc.frac_df=2.3D0;

% Physicochemical PropertiprFloc.es of Solid and Liquid
prFloc.par_den=1.80D3;
prFloc.wat_den=1.05D3;
prFloc.g=9.81D0;
prFloc.vmu=1.002D-3;
prFloc.bolz=1.38D-23;
prFloc.temp=2.93D2;

% Constants (make 'MKS' units)- Check Spicer & Pratsinis (AICHE, 1996)
prFloc.pdia=15.0D-6;
prFloc.ffdia=100.0D-6;
prFloc.cgel=100.0d0;
prFloc.pri_vol=1.0D0/6.0D0*3.1416d0*prFloc.pdia^3;

% Initial Concentration
prIni.pconc=0.32D0;
prIni.pnconc=prIni.pconc/prFloc.par_den/prFloc.pri_vol;

% Erosion/Deposition Constants
prBot.ers_m=0.25D-3;
prBot.tau_c=1.0D0;
prBot.fntc=1.0D0/6.0D0*3.1416D0*prFloc.pdia^3*prFloc.par_den;

% coefficient for boundary conditions for velocity and turbulent equations
prFlow.v_karm=0.4d0;

% Grid Information
prGrid.t_height = 7.0d0;
prGrid.dz = 1.0d0;
prGrid.np = int32(prGrid.t_height/prGrid.dz);

% time step size
prTime.dt=10.0d0;
printHour = 0.5; %hours
prTime.print= printHour*3600/prTime.dt;

% Now save in plugin structure
plugin.prFloc = prFloc;
plugin.prIni = prIni;
plugin.prBot = prBot;
plugin.prFlow = prFlow;
plugin.prGrid = prGrid;
plugin.prTime = prTime;

% Now initialize beta for explicit solution
plugin.beta.pbe1 = zeros(prGrid.np,1); 
plugin.beta.bm1 = zeros(prGrid.np,1);
plugin.beta.sh1 = zeros(prGrid.np,1); 
plugin.beta.pbe2 = zeros(prGrid.np,1);
plugin.beta.bm2 = zeros(prGrid.np,1); 
plugin.beta.ds2 = zeros(prGrid.np,1);
plugin.beta.sh2 = zeros(prGrid.np,1); 
plugin.beta.pbe3 = zeros(prGrid.np,1);
plugin.beta.bm3 = zeros(prGrid.np,1); 
plugin.beta.ds3 = zeros(prGrid.np,1);
plugin.beta.sh3 = zeros(prGrid.np,1);
% Ininitialize ae2 and ae3 for implicit solution
plugin.ae2 = zeros(prGrid.np,1); 
plugin.ae3 = zeros(prGrid.np,1);

plugin.delz = ones(prGrid.np,1) * ...   % delz of settling column
    prGrid.dz;   
plugin.np = prGrid.np;          	% # nodes
plugin.ts_onemab = zeros(1,5);          % initial condition (intial seed)
plugin.sig_b = 0.0d0;                   % initial sig_b
plugin.time = 0.0d0;                    % initial time

end
