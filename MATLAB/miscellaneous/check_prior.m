function r_arg = check_prior(Par_info,DREAMPar,M)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Determines whether the prior handle returns a pdf or logpdf           %%
%%                                                                       %%
%% SYNOPSIS: r_arg = check_prior(Par_info,DREAMPar)                      %%
%%  where                                                                %%
%%   Par_info  [input] Parameter structure                               %%
%%   DREAMPar  [input] DREAM Package structure                           %%
%%   M         [input] OPT: Number of samples                            %%
%%   r_arg     [outpt] string: 'pdf' or 'logpdf'                         %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, July 2024                               %%
%% University of California Irvine                                       %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3, M = 100; end

r_arg = 'pdf';               % Default is pdf

% Draw samples from prior and then evaluate prior handle of pdf (or logpdf) 
% Return argument helps decide whether it is a pdf or logpdf 
% This may not work for highly peaked priors as pdf > 1 too often
if strcmp(Par_info.u,'yes')
    PrX = nan(M,DREAMPar.d);
    for zz = 1:M
        for j = 1:DREAMPar.d
            PrX(zz,j) = Par_info.prior{j}(Par_info.prior_rnd{j}(1));
        end
    end
%    if sum(sum(PrX>=0)) > M*DREAMPar.d/2
    if sum(sum(PrX<0)) > 0
        r_arg = 'logpdf';
    end
else
    PrX = nan(M,1);
    for zz = 1:M
        PrX(zz,1) = Par_info.prior(Par_info.prior_rnd(1));
    end
%    if sum(PrX>=0) > M/2
    if sum(PrX<0) > 0
        r_arg = 'logpdf';
    end
end

end