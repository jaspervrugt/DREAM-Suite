function [FX,S] = Evaluate_target(X,DREAMPar,Meas_info,options, ...
    f_handle,verbose)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function evaluates the function handle f_handle                  %%
%%                                                                       %%
%% SYNOPSIS: [FX,S] = Evaluate_target(X,DREAMPar,Meas_info,options,...   %%
%%                        f_handle,verbose)                              %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, June 2006                               %%
%% Los Alamos National Laboratory                                        %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 6
    verbose = 0;
else
    clear('prt_progress');      % Waitbar initialization
end

N = size(X,1);                  % # parameter vectors?

% Preallocate FX
if strcmp(options.modout,'yes') && Meas_info.n > 0
    switch strcmp(options.DB,'yes')                 
        case 1
            FX = nan(Meas_info.n+Meas_info.n_S,N);
        otherwise
            FX = nan(Meas_info.n,N);
    end
end

% Evaluate target (= model)
switch DREAMPar.CPU                                 
    case 1    % Sequential evaluation
        if isfield(DREAMPar,'case')                 % Forward model vectorized/not
            FX = f_handle(X);                       % Pass x directly to function
        else
            for ii = 1:N                            % Otherwise, candidate points separately
                FX(:,ii) = f_handle(X(ii,:));
                if verbose                          % Print output
                    prt_progress(DREAMPar,N);
                end
            end
        end
    otherwise   % Parallel evaluation
        if strcmp(options.IO,'yes')
            EXAMPLE_dir = pwd;
            parfor ii = 1:N                         % Loop over candidate points/workers
                task = getCurrentTask();            % Determine task
                id = get(task, 'ID');               % Determine work ID
                cd([EXAMPLE_dir,'/',num2str(id)]);  % Right directory 
                FX(:,ii) = f_handle( X(ii,:) );     %#ok Evaluate target/model
                if verbose && id == 1               % Print output
                    prt_progress(DREAMPar,N);
                end
            end
            cd(EXAMPLE_dir)
        elseif strcmp(options.IO,'no')
            parfor ii = 1:N                         % Loop over candidate points/workers
                FX(:,ii) = f_handle(X(ii,:));       %#ok evaluate target/model
                if verbose                          % Print output
                    task = getCurrentTask();        % Determine task
                    id = get(task, 'ID');           % Determine work ID
                    if id == 1                      % Print output
                        prt_progress(DREAMPar,N); 
                    end                             
                end
            end
        end
end

% Now check if diagnostic Bayes is used
switch strcmp(options.DB,'yes')
    case 1
        S = FX(Meas_info.n + 1 : end , 1 : N);      % Extract summary metrics
        FX = FX(1 : Meas_info.n , 1 : N);           % Remove from fx values
    otherwise
        S = [];
end

% Go to next line if verbose is activated
if verbose, fprintf('\n'); fprintf('Posterior simulation ... done\n'); end

end
