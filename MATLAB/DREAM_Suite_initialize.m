function [chain,output,X,FX,SX,Table_gamma,CR,pCR,tCR,TsdX_Xp,sdX_Xp,...
    loglik,EX,std_mX,id_r,id_c,iloc,iter,gen,LLX,MAP_info] = ...
    DREAM_Suite_initialize(method,DREAMPar,Par_info,Meas_info,...
    Lik_info,options,f_handle,MAP_info)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%             DDDDD     RRRRR     EEEEEEE     AAA     MM    MM            %
%             DDDDDD    RRRRRR    EEEEEEE    AAAAA    MM    MM            %
%             DD  DD    RR   RR   EE        AA   AA   MMM  MMM            %
%             DD   DD   RR  RR    EEEE      AA   AA   MMMMMMMM            %
%             DD   DD   RRRRR     EEEE      AAAAAAA   MMM  MMM            %
%             DD  DD    RR RR     EE        AAAAAAA   MM    MM            %
%             DDDDDD    RR  RR    EEEEEEE   AA   AA   MM    MM            %
%             DDDDD     RR   RR   EEEEEEE   AA   AA   MM    MM            %
%                                                                         %
%              SSSSSSSS  UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SSSSSSS   UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SS        UU    UU   II       TT       EE                  %
%              SSSS      UU    UU   II       TT       EEEE                %
%                 SSSS   UU    UU   II       TT       EEEE                %
%                   SS   UU    UU   II       TT       EE                  %
%               SSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%              SSSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Initializes main variables and chains of DREAM-Suite                    %
%                                                                         %
% SYNOPSIS: [chain,output,X,FX,S,Table_gamma,CR,pCR,tCR,TsdX_Xp, ...      %
%               sdX_Xp,loglik,EX,std_mX,iloc,iter,gen,LLX,MAP_info] = ... %
%               DREAM_Suite_initialize(method,DREAMPar,Par_info, ...      %
%               Meas_info,Lik_info,options,f_handle,MAP_info)             %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

switch method
    case {'dream_zs','dream_dzs','dream_kzs','mtdream_zs'}
        N = DREAMPar.m0 + DREAMPar.N;
    case {'dream','dream_d'}
        N = DREAMPar.N;
end
% Initialize X
X = nan(N,DREAMPar.d);

% Generate the initial states (and/or initial archive) of the chains
switch Par_info.initial
    case {'uniform'}    % Initialize chains from multivariate uniform
        X = repmat(Par_info.min,N,1) + rand(N,DREAMPar.d) .* ...
            ( repmat(Par_info.max - Par_info.min,N,1) );
    case {'latin'}      % Initialize chains with Latin hypercube sampling
        X = LH_sampling(Par_info.min,Par_info.max,N);
    case {'normal'}     % Initialize chains with (multi)-normal distribution
        X = repmat(Par_info.mu,N,1) + randn(N,DREAMPar.d) * chol(Par_info.cov);
    case {'prior'}      % Initialize chains by drawing from prior distribution
        switch Par_info.u
            case {'yes'}
                % Univariate prior: Draw one parameter at a time
                for qq = 1:DREAMPar.d
                    for zz = 1:N
                        X(zz,qq) = Par_info.prior_rnd{qq}(1);
                    end
                end
            case {'no'}
                % Multivariate prior: Draw all parameters at once
                for zz = 1:N
                    X(zz,1:DREAMPar.d) = Par_info.prior_rnd(1);
                end
        end
    case {'user'}       % Initial state of chains specified by user
        for zz = 1:N
            X(zz,1:DREAMPar.d) = Par_info.x0(zz,1:DREAMPar.d);
        end
end

% If normalized sampling: ranges of X between 0 and 1
switch Par_info.norm
    case 1 % X must be between 0 and 1
        if isfield(Par_info,'minun')
            X = bsxfun(@rdivide,X - Par_info.minun,Par_info.maxun - ...
                Par_info.minun);
        else
            X = bsxfun(@rdivide,X - min(X),max(X)-min(X));
            % --> this may cause problems as 0 - 1 bounds are determined by
            %     prior draws
        end
    otherwise % Everything OK
end

% If specified do boundary handling ( "Bound","Reflect","Fold")
if isfield(Par_info,'boundhandling')
    [X,v] = Boundary_handling(X,Par_info);
    % JAV, June 12, 2018: Make sure to verify "v" -> violate vector; important for
    % initial X and the starting samples of Z (not yet checked/implemented here)
else
    v = ones(N,1) < 0;
end

if any(strcmp(method,{'dream_d','dream_dzs'}))
    % Now transform X to discrete space
    X = Discrete_space(X,Par_info);
end

% Now check method to use
switch method
    case {'dream_zs','dream_dzs','dream_kzs','mtdream_zs'}
        % Define initial archive Z; prior/likelihood not required
        Z = [ X(1:DREAMPar.m0,1:DREAMPar.d) NaN(DREAMPar.m0,2) ];
        % Now open binary file
        fid_Z = fopen('Z.bin','w+','n');
        % Now write Z' to binary file; Z' is required for fseek later
        fwrite(fid_Z,Z','double');
        % Now close file
        fclose(fid_Z);
        % Now remove from X the Z samples to have left initial population
        X = X(DREAMPar.m0 + 1:DREAMPar.m0 + DREAMPar.N,1:DREAMPar.d);
        % Do the same thing with vector v
        v = v(DREAMPar.m0 + 1:DREAMPar.m0 + DREAMPar.N,1);
        % And delete the content of the FXZ.bin file
        fid_FXZ = fopen('FXZ.bin','w'); fclose(fid_FXZ);
        % To perfectly synchronize Z.bin and FXZ.bin one should fill
        % FXZ.bin with nan for the first m0 rows
    case {'dream','dream_d'}
        % do nothing - no used of past samples
end

if strcmp(method,'dream_kzs')
    % Now evaluate the model ( = pdf ) and return FX
    [FXZ] = Evaluate_target(Z(:,1:DREAMPar.d),DREAMPar,Meas_info,options,f_handle);
    % Now save output
    fid_FXZ = fopen('FXZ.bin','w+','n');
    % Now write to file
    fwrite(fid_FXZ,FXZ,'double');
    % Now close file
    fclose(fid_FXZ);
end

%% Compute sandwich properties - for a(theta) correction
if isfield(MAP_info,'map')
    map_un = X_unnormalize(MAP_info.map,Par_info);
    % Now evaluate the model ( = pdf ) and return FX
    [Fmap,~] = Evaluate_target(map_un,DREAMPar,Meas_info,...
        options,f_handle);
    % Compute maximum of log-likelihood - without sandwich adjustment
    [~,~,~,~,~,LLmax] = Calc_likelihood(map_un,Fmap,DREAMPar,...
        Par_info,Meas_info,Lik_info,options);
    % Store data of MAP solution, likelihood maximum, parameter values
    MAP_info.loglik = sum(LLmax); MAP_info.map_un = map_un;
    % Fisher and Godambe information
    %       MAP_info.Fisher = MAP_info.An;
    %       MAP_info.Godambe = MAP_info.An * inv(MAP_info.Betan) * MAP_info.An;
    MAP_info.Fisher = MAP_info.An; MAP_info.Godambe = (MAP_info.An / ...
        MAP_info.Betan) * MAP_info.An;
    % Display MAP_info on screen
    disp(MAP_info)
end

% Now normalize parameter values or not?
X_un = X_unnormalize(X,Par_info);
% Now evaluate the model ( = pdf ) and return FX
[FX,SX] = Evaluate_target(X_un,DREAMPar,Meas_info,options,f_handle);
% NEW: Compute the log(prior) of candidate points
logPR_X = Eval_prior(X_un,SX,Par_info,Meas_info,options);
% NEW: Compute the log-likelihood of candidate points
[logL_X,EX,std_mX,~,~,LLX] = Calc_likelihood(X_un,FX,DREAMPar,...
    Par_info,Meas_info,Lik_info,options,MAP_info);
% NEW: Now append log(prior) and log-likelihood to matrix Xp
X(1:DREAMPar.N,DREAMPar.d+1:DREAMPar.d+2) = [logPR_X logL_X];
% prop. loglik if range violation: Par_info.boundhandling = 'reject'
% then, Xp(v==1,:) [= Xp(v,:)] will be declined
X(v,DREAMPar.d+2) = -inf;

% Store the model simulations (if appropriate)
DREAM_Suite_store_results ( options , [ FX ; SX]  , Meas_info , 'w+' );

% Initialize matrix with log_likelihood of each chain
loglik = nan(DREAMPar.T,DREAMPar.N+1);
% Save history log density of individual chains
loglik(1,1:DREAMPar.N+1) = [ DREAMPar.N X(1:DREAMPar.N,DREAMPar.d+2)' ];

% Define number of elements of several fields of structure output
nr = floor(DREAMPar.T/DREAMPar.steps);
% Define structure output
output = struct('AR',nan(nr,2),'MR_stat',nan(nr,2),'R_stat',...
    nan(nr,DREAMPar.d+1),'CR',nan(nr,DREAMPar.nCR+1),'outlier',[]);
% Now set first line of acceptance rate
output.AR(1,1) = DREAMPar.N;

% Define selection probability of each crossover
pCR = (1/DREAMPar.nCR) * ones(1,DREAMPar.nCR);
% Sample randomly the crossover values [ 1 : DREAMPar.nCR ]
CR = reshape(randsample(1:DREAMPar.nCR,DREAMPar.N * DREAMPar.steps,...
    'true',pCR),DREAMPar.N,DREAMPar.steps);
% Initialize the tCR values, TsdX_Xp and sdX_Xp values
tCR = ones(1,DREAMPar.nCR); TsdX_Xp = ones(1,DREAMPar.nCR);
sdX_Xp = zeros(DREAMPar.N,DREAMPar.steps);
% Save pCR values in memory
output.CR(1,1:DREAMPar.nCR+1) = [ DREAMPar.N pCR ];

% Initialize array (3D-matrix) of chain trajectories
chain = nan(floor(DREAMPar.T/DREAMPar.thinning),DREAMPar.d+2,DREAMPar.N);
% Set initial state of chains equal to the first sampled X values
chain(1,1:DREAMPar.d+2,1:DREAMPar.N) = ...
    reshape(X',1,DREAMPar.d+2,DREAMPar.N);

% Generate Table with jump rates: more efficient to read from Table
Table_gamma = nan(DREAMPar.d,DREAMPar.delta);
for zz = 1:DREAMPar.delta
    Table_gamma(1:DREAMPar.d,zz) = 2.38./sqrt(2 * zz * (1:DREAMPar.d)');
end
% Reduce/increase linearly jumprate if so desired by user
Table_gamma = DREAMPar.beta0 * Table_gamma;

% Initialize few important counters
iloc = 1; iter = 2; gen = 2;
% How many total candidate points in all chains
Nmt = DREAMPar.mt * DREAMPar.N;
% Index of reference and selected/center point
id_r = 1 : Nmt; id_c = DREAMPar.mt : DREAMPar.mt : Nmt;
% Remove center point from idx_r
id_r(id_c) = [];

% Print to screen the choice of measurement error function
if isfield(Meas_info,'sigma2') || ~isempty(Meas_info.Sigma)
    switch Lik_info.method
        case 0
            fprintf('DREAM_Suite: s0/s1 not considered \n');
        case 1
            fprintf(['DREAM_Suite: CONSTANT measurement error: s0 ' ...
                'fixed at default value\n']);
        case 2
            fprintf(['DREAM_Suite: CONSTANT measurement error: s0 ' ...
                'treated as phantom variable\n']);
        case 3
            fprintf(['DREAM_Suite: CONSTANT measurement error: s0 ' ...
                'estimated\n']);
        case 4
            fprintf(['DREAM_Suite: CONSTANT measurement error: s0 ' ...
                'estimated\n']);
        case 5
            fprintf(['DREAM_Suite: NONCONSTANT measurement error: s0' ...
                ' fixed at default value and s1 treated as phantom ' ...
                'variable\n']);
        case 6
            fprintf(['DREAM_Suite: NONCONSTANT measurement error: s0' ...
                ' fixed at default value and s1 estimated\n']);
        case 7
            fprintf(['DREAM_Suite: NONCONSTANT measurement error: s0' ...
                ' estimated and s1 treated as phantom variable\n']);
        case 8
            fprintf(['DREAM_Suite: NONCONSTANT measurement error: s0' ...
                ' estimated and s1 estimated\n']);
    end
end

% % fprintf('\n');
% % fprintf(['  ---------------------------------------------------------' ...
% %     '--------------------------------------------------  \n']);
% % fprintf('  DDDDD    RRRRRR   EEEEEE     A     M     M    PPPP        A     CCCCCCC  K      KK     A     GGGGGG  EEEEEE \n');
% % fprintf('  DDDDDD   RRRRRRR  EEEEEE    AAA    M     M    PPPPP      AAA    CCCCCC   K     KK     AAA    G    G  EEEEEE \n');
% % fprintf('  DDDDDDD  R     R  E        AAAAA   MM   MM    PPPPPP    AAAAA   CC       K     KK    AAAAA   G    G  E      \n');
% % fprintf('  D     D  R    R   E       A     A  MMM MMM    P    PP  A     A  C        K    KK    A     A  G    G  E      \n');
% % fprintf('  D     D  R   R    EEE     A     A  MMMMMMM    P    PP  A     A  C        KKKKKK     A     A  GGGGGG  EEE    \n');
% % fprintf('  D     D  RRRRR    EEE     AAAAAAA  M MMM M    PPPPPP   AAAAAAA  C        KKKKKK     AAAAAAA  GGGGGG  EEE    \n');
% % fprintf('  D     D  RRRRR    E       AAAAAAA  M  M  M    P        AAAAAAA  C        K    KK    AAAAAAA       G  E      \n');
% % fprintf('  DDDDDDD  R   R    E       A     A  M     M    P        A     A  CC       K     KK   A     A       G  E      \n');
% % fprintf('  DDDDDD   R    R   EEEEEE  A     A  M     M    P        A     A  CCCCCC   K     KK   A     A      GG  EEEEEE \n');
% % fprintf('  DDDDD    R     R  EEEEEE  A     A  M     M    P        A     A  CCCCCCC  K      KK  A     A  GGGGGG  EEEEEE \n');
% % fprintf(['  ---------------------------------------------------------' ...
% %     '--------------------------------------------------  \n']);
% % fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
% % fprintf('\n');

% Upper case
fprintf('\n');
fprintf(['  ---------------------------------------------------------' ...
    '--------------------------------------------------  \n']);
fprintf('  DDDDD    RRRRRR   EEEEEE     A     M     M    PPPPPP      A     CCCCCCC  K      K     A      GGGGGG  EEEEEE \n');
fprintf('  D    D   R     R  E         A A    M     M    P     P    A A    C        K     K     A A    G     G  E      \n');
fprintf('  D     D  R     R  E        A   A   M     M    P     P   A   A   C        K     K    A   A   G     G  E      \n');
fprintf('  D     D  R    R   E       A     A  MM   MM    P     P  A     A  C        K    K    A     A  G     G  E      \n');
fprintf('  D     D  R   R    EEEEE   A     A  M M M M -- PPPPPP   A     A  C        KKKKK     A     A   GGGGGG  EEEEE  \n');
fprintf('  D     D  RRRRR    E       AAAAAAA  M  M  M -- P        AAAAAAA  C        K    K    AAAAAAA        G  E      \n');
fprintf('  D     D  R    R   E       A     A  M     M    P        A     A  C        K     K   A     A        G  E      \n');
fprintf('  D    D   R     R  E       A     A  M     M    P        A     A  C        K     K   A     A       GG  E      \n');
fprintf('  DDDDD    R     R  EEEEEE  A     A  M     M    p        A     A  CCCCCCC  K      K  A     A   GGGGGG  EEEEEE \n');
fprintf(['  ---------------------------------------------------------' ...
    '--------------------------------------------------  \n']);
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n');

% fprintf('\n');
% fprintf(['  ---------------------------------------------------------' ...
%     '--------------------------------------------------  \n']);
% fprintf(['  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><' ...
%     '><><><><><><><><><><><><><><><><><><><><><><><><><  \n']);
% 
switch method
    case 'dream'
        fprintf('%37s ddddd   rrrrr  eeeeee   aa   m    m         \n','');
        fprintf('%37s d    d  r    r e       a  a  m    m         \n','');
        fprintf('%37s d     d r    r e      a    a mm  mm         \n','');
        fprintf('%37s d     d r    r eeeee  a    a mmmmmm         \n','');
        fprintf('%37s d     d rrrrrr e      aaaaaa mm  mm         \n','');
        fprintf('%37s d     d r    r e      a    a m    m         \n','');
        fprintf('%37s dddddd  r    r eeeeee a    a m    m         \n','');
    case 'dream_d'
        fprintf('%34s ddddd   rrrrr  eeeeee   aa   m    m            \n','');
        fprintf('%34s d    d  r    r e       a  a  m    m            \n','');
        fprintf('%34s d     d r    r e      a    a mm  mm            \n','');
        fprintf('%34s d     d r    r eeeee  a    a mmmmmm            \n','');
        fprintf('%34s d     d rrrrrr e      aaaaaa mm  mm ddddd      \n','');
        fprintf('%34s d     d r    r e      a    a m    m d    d     \n','');
        fprintf('%34s dddddd  r    r eeeeee a    a m    m d    d     \n','');
	    fprintf('%34s                                     d    d     \n','');
	    fprintf('%34s                                     ddddd	     \n','');
    case 'dream_zs'
        fprintf('%30s ddddd   rrrrr  eeeeee   aa   m    m               \n','');
        fprintf('%30s d    d  r    r e       a  a  m    m               \n','');
        fprintf('%30s d     d r    r e      a    a mm  mm               \n','');
        fprintf('%30s d     d r    r eeeee  a    a mmmmmm               \n','');
        fprintf('%30s d     d rrrrrr e      aaaaaa mm  mm zzzzzz ssssss \n','');
        fprintf('%30s d     d r    r e      a    a m    m    zz  s      \n','');
        fprintf('%30s dddddd  r    r eeeeee a    a m    m   zz   ssssss \n','');
	    fprintf('%30s                                      zz         s \n','');
	    fprintf('%30s                                     zzzzzz ssssss \n','');
    case 'dream_dzs'
   	    fprintf('%26s ddddd   rrrrr  eeeeee   aa   m    m                      \n','');
        fprintf('%26s d    d  r    r e       a  a  m    m                      \n','');
        fprintf('%26s d     d r    r e      a    a mm  mm                      \n','');
        fprintf('%26s d     d r    r eeeee  a    a mmmmmm                      \n','');
        fprintf('%26s d     d rrrrrr e      aaaaaa mm  mm ddddd  zzzzzz ssssss \n','');
        fprintf('%26s d     d r    r e      a    a m    m d    d    zz  s      \n','');
        fprintf('%26s dddddd  r    r eeeeee a    a m    m d    d   zz   ssssss \n','');
	    fprintf('%26s                                     d    d  zz         s \n','');
	    fprintf('%26s                                     ddddd  zzzzzz ssssss \n','');
    case 'mtdream_zs'
        fprintf('%21s m    m ttttttt    ddddd   rrrrr  eeeeee   aa   m    m               \n','');
        fprintf('%21s m    m    t       d    d  r    r e       a  a  m    m               \n','');
        fprintf('%21s mm  mm    t       d     d r    r e      a    a mm  mm               \n','');
        fprintf('%21s mmmmmm    t    -- d     d r    r eeeee  a    a mmmmmm               \n','');
        fprintf('%21s mm  mm    t       d     d rrrrrr e      aaaaaa mm  mm zzzzzz ssssss \n','');
        fprintf('%21s m    m    t       d     d r    r e      a    a m    m    zz  s      \n','');
        fprintf('%21s m    m    t       dddddd  r    r eeeeee a    a m    m   zz   ssssss \n','');
	    fprintf('%21s                                                        zz         s \n','');
	    fprintf('%21s                                                       zzzzzz ssssss \n','');
end

% fprintf(['  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><' ...
%     '><><><><><><><><><><><><><><><><><><><><><><><><><  \n']);
% fprintf('\n');

end