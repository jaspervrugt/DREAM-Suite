function [TsdX_Xp,cCR,pCR] = Calc_pCR(DREAMPar,sdX_Xp,TsdX_Xp,cCR,CR)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function updates the selection probabilities of the nCR crossover values      %%
%%                                                                                    %%
%% SYNOPSIS: [TsdX_Xp,cCR,pCR] = Calc_pCR(DREAMPar,sdX_Xp,TsdX_Xp,cCR,CR)             %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, June 2006                                            %%
%% Los Alamos National Laboratory                                                     %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

for i = 1:DREAMPar.N
    for j = 1:DREAMPar.steps
        % Traveled standardized Euclidean distance due to CR(i,j)
        TsdX_Xp(CR(i,j)) = TsdX_Xp(CR(i,j)) + sdX_Xp(i,j);
        % Count of CR(i,j) [= how many times has this value been used]
        cCR(CR(i,j)) = cCR(CR(i,j)) + 1;
    end
end
pCR = TsdX_Xp ./ cCR;   % New selection probability crossover values
pCR = pCR ./ sum(pCR);  % Normalize pCR so that sum(pCR) is equal to one

end