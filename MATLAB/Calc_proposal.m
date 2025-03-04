function [Xp,log_alfa_sn,CR,v] = Calc_proposal(method,X,EX,std_mX,CR,...
    DREAMPar,Table_gamma,Par_info,Meas_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function generates candidate points using a discrete proposal    %%
%%  distribution                                                         %%
%%                                                                       %%
%% SYNOPSIS: [Xp,log_alfa_sn,CR,v] = Calc_proposal(method,X,EX,...       %%
%%              std_mX,CR,DREAMPar,Table_gamma,Par_info,Meas_info)       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Feb 2007                                %%
%% Los Alamos National Laboratory                                        %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Set N jump vectors and candidate points equal to zero
[Xp,dX] = deal(zeros(DREAMPar.N,DREAMPar.d));
% Define log of alfa snooker
log_alfa_sn = zeros(DREAMPar.N,1);
% Determine to do parallel direction or snooker update
jump_method = randsample(1:3,1,'true',[DREAMPar.pparallel ...
    DREAMPar.psnooker DREAMPar.pkalman]);

switch jump_method
    case {1,2}
        switch method
            case {'dream','dream_d'} % Current chain states for proposals
                % Randomly permute numbers [1,...,N-1] N times
                [~,draw] = sort(rand(DREAMPar.N-1,DREAMPar.N));
                % Set current X values equal to Z
                Z = X;
            case {'dream_zs','dream_dzs','dream_kzs'}
                % Permute [1,...,N-1] N times so that Z samples differ
                draw = repmat(1:6,DREAMPar.N,1)';
                % Without replacement draw rows from Z to create proposals
                R = randsample(DREAMPar.m,DREAMPar.select); 
                Z = nan(DREAMPar.select,DREAMPar.d);
                % Open the external archives with historical chain states
                fid_Z = fopen('Z.bin');
                % Now read from file
                for zz = 1 : DREAMPar.select
                    % Open file "Z.bin"
                    fseek(fid_Z,(R(zz)-1) * (DREAMPar.d + 2) * 8,'bof');
                    % Now load row and store in matrix Z
                    Z(zz,1:DREAMPar.d) = fread(fid_Z,DREAMPar.d,'double')';
                end
                % Close the file
                fclose(fid_Z);
        end
    case {3} % dream_kzs: Past states and their simulations for proposals
        % Without replacement draw rows from Z to create proposals
        R = randsample(DREAMPar.m,DREAMPar.M); Z = nan(DREAMPar.M,DREAMPar.d);
        % Open external archives of past chain states and their simulations
        fid_Z = fopen('Z.bin'); fid_FXZ = fopen('FXZ.bin');
        % Initialize FXZ
        FXZ = nan(DREAMPar.M,Meas_info.n);
        % Now read from file
        for zz = 1 : DREAMPar.M
            % Now seek the right rows of "Z.bin"
            fseek(fid_Z,(R(zz)-1) * (DREAMPar.d + 2) * 8,'bof');
            % Now load row and store in matrix Z
            Z(zz,1:DREAMPar.d) = fread(fid_Z,DREAMPar.d,'double')';
            % Now seek the right rows of "FXZ.bin"
            fseek(fid_FXZ,(R(zz)-1) * Meas_info.n * 8,'bof');
            % Now load row and store in matrix FXZ
            FXZ(zz,1:Meas_info.n) = fread(fid_FXZ,Meas_info.n,'double')';
        end
        % Close the two files
        fclose(fid_Z); fclose(fid_FXZ);
end

% Determine how many chain pairs to use for each individual chain
delta = randsample( 1 : DREAMPar.delta , DREAMPar.N , 'true' , ...
    (1/DREAMPar.delta) * ones(1,DREAMPar.delta) )';
% Uniform random numbers each chain to determine which dimensions to update
rnd_cr = rand(DREAMPar.N,DREAMPar.d);
% Ergodicity perturbation
eps = DREAMPar.zeta * randn(DREAMPar.N,DREAMPar.d);
% Ergodicity for each individual chain
rnd_jump = DREAMPar.lambda * (2 * rand(DREAMPar.N,DREAMPar.d) - 1); 
%rnd_jump = DREAMPar.lambda * (2 * rand(DREAMPar.N,1) - 1);

switch jump_method
    case 1  % PARALLEL DIRECTION PROPOSAL DISTRIBUTION
        % Determine when jumprate is 1
        gamma = randsample([0 1],DREAMPar.N,true,...
            [ 1 - DREAMPar.p_unit_gamma DREAMPar.p_unit_gamma ]);
        % Create N proposals
        for i = 1:DREAMPar.N
            % Derive vectors a and b
            a = DREAMPar.R(i,draw(1:delta(i),i));
            b = DREAMPar.R(i,draw(delta(i)+1:2*delta(i),i));
            % Which gamma to use?
            if ( gamma(i) == 1 )
                % If gamma = 1 --> full dimensional update (mode jumping)
                A = 1:DREAMPar.d; gamma_d = 1; CR(i,1) = DREAMPar.nCR;
            else
                % Derive subset A with dimensions to sample
                A = find( rnd_cr(i,1:DREAMPar.d) < CR(i,1)/DREAMPar.nCR );
                % How many dimensions are sampled?
                d_prime = numel(A);
                % Make sure that at least one dimension is selected!
                if ( d_prime == 0 )
                    idx = randperm(DREAMPar.d); A = idx(1); d_prime = 1;
                end
                % Unpack jump rate
                gamma_d = Table_gamma(d_prime,delta(i)); 
                % gamma_d = (0.5 + rand/2) * sqrt(1/DE_pairs(i));
            end
            % Calculate jump
            dX(i,A) = ( 1 + rnd_jump(i,A) ) .* gamma_d .* ...
                sum( Z(a,A) - Z(b,A) , 1 ) + eps(i,A);
            % Compute candidate point
            Xp(i,1:DREAMPar.d) = X(i,1:DREAMPar.d) + dX(i,1:DREAMPar.d);
        end
    case 2  % SNOOKER PROPOSAL DISTRIBUTION
        % Loop over the individual chains
        for i = 1:DREAMPar.N
            % Sample a, b and c from external archive Z
            a = DREAMPar.R(i,1); b = DREAMPar.R(i,2); c = DREAMPar.R(i,3);
            % Projection vector X(i,:) - Za
            F = X(i,1:DREAMPar.d) - Z(c,1:DREAMPar.d); 
            % Make sure inner product F*F' is at least realmin 
            FF = max(F*F',realmin);
            % Orthogonally projection of zR1 and zR2 onto F
            zp = F * (sum((Z(a,1:DREAMPar.d) - Z(b,1:DREAMPar.d)).*F)/FF);
            % Define JumpRate -- uniform rand number between 1.2 and 2.2
            gamma_s = 1.2 + rand;
            % And define the jump
            dX(i,1:DREAMPar.d) = ( 1 + rnd_jump(i,1:DREAMPar.d) ) .* ...
                gamma_s .* zp + eps(i,1:DREAMPar.d);
            % Compute candidate points inside loop for snooker correction
            Xp(i,1:DREAMPar.d) = X(i,1:DREAMPar.d) + dX(i,1:DREAMPar.d);
            % First calculate Euclidean distance of Xp to Zc ( = L2 norm )
            XpZ = norm(Xp(i,1:DREAMPar.d) - Z(c,1:DREAMPar.d),2) + realmin;
            % Then calculate Euclidean distance of X to Zc ( = L2 norm )
            XZ = norm(X(i,1:DREAMPar.d) - Z(c,1:DREAMPar.d),2) + realmin;
            % Now calculate snooker correction (in log space)
            log_alfa_sn(i,1) = (DREAMPar.d - 1) * log(XpZ/XZ);
            % norm(x-z,2) = Euclid. dist.: ||.|| --> sqrt( (x1-z1)^2 + 
            %     (x2-z2)^2 + ... + (xn-zn)^2 )
            % alfa_s = ( || Xp - Zc || / || X - Zc || ).^((DREAMPar.d - 1))
        end
    case 3  % KALMAN PROPOSAL DISTRIBUTION
        % Replicate mean of sampled vectors of Z
        m_Z = repmat(mean(Z,1),DREAMPar.M,1);
        % Replicate mean of model output of sampled vectors of Z
        m_FX = repmat(mean(FXZ,1),DREAMPar.M,1);
        % Compute cross-covariance between Z and FXZ
        C_ZFXZ = (Z - m_Z)' * (FXZ - m_FX)/(DREAMPar.M - 1);
        % Compute auto-covariance of model output
        C_FXZFXZ = (FXZ - m_FX)' * (FXZ - m_FX)/(DREAMPar.M - 1);
        for i = 1:DREAMPar.N
            if ~isfield(Meas_info,'R')
                if all(all(isnan(std_mX), 1))
                    R = 0;
                    % Kalman gain - without measurement error
                    % K = C_ZFXZ/C_FXZFXZ;
                    % Vector r is a zero-vector
                    r = zeros(Meas_info.n,1);
                else
                    R = diag(std_mX(1:Meas_info.n,i).^2);
                    % Kalman gain: measurement error from Sigma
                    % K = C_ZFXZ/(C_FXZFXZ + cov(std_mX(1:Meas_info.n,i)) ...
                    %     * eye(Meas_info.n));
                    % Vector r from measurement error model
                    r = randn(Meas_info.n,1) .* std_mX(1:Meas_info.n,i);
                end
            elseif ~isempty(Meas_info.R)
                % Compute Kalman gain (R known)
                R = Meas_info.R;
                % K = C_ZFXZ/(C_FXZFXZ + Meas_info.R);
                if all(all(isnan(std_mX), 1))
                    % Vector r from measurement error covariance matrix
                    r = mvnrnd(zeros(1,Meas_info.n),R)';
                else
                    % Vector r from measurement error model
                    r = randn(Meas_info.n,1) .* std_mX(1:Meas_info.n,i);
                end
            else % Meas_info.R = []
                if all(all(isnan(std_mX), 1))
                    R = 0;
                    % Compute Kalman gain: without measurement error
                    % K = C_ZFXZ/C_FXZFXZ;
                    % Vector r is a zero-vector
                    r = zeros(Meas_info.n,1);
                else
                    R = diag(std_mX(1:Meas_info.n,i).^2);
                    % Compute Kalman gain: measurement error from Sigma
                    % K = C_ZFXZ/(C_FXZFXZ + cov(std_mX(1:Meas_info.n,i)) ...
                    %     * eye(Meas_info.n));
                    % Vector r from measurement error model
                    r = randn(Meas_info.n,1) .* std_mX(1:Meas_info.n,i);
                end
            end
            % Kalman gain: 
            K = C_ZFXZ / (C_FXZFXZ + R);
            % Kalman jump
            dx_K = K * ( r  - EX(1:Meas_info.n,i) );
            % Jump vector
            dX(i,:) = ( 1 + rnd_jump(i,1:DREAMPar.d) ) .* dx_K' + ...
                eps(i,1:DREAMPar.d);
            % Candidate point
            Xp(i,1:DREAMPar.d) = X(i,1:DREAMPar.d) + dX(i,1:DREAMPar.d);
        end
end

% Boundary handling ('bound'/'reflect'/'fold'/'reject')
if isfield(Par_info,'boundhandling')
    [Xp,v] = Boundary_handling(Xp,Par_info);
else
    % Define N x 1 violation vector: logical 0) in bound; 1) out of bound
    v = zeros(DREAMPar.N,1) > 1;
end
if any(strcmp(method,{'dream_d','dream_dzs'}))
    % Now transform Xp to discrete space
    Xp = Discrete_space(Xp,Par_info);
end

end
