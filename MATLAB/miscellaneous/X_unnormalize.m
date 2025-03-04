function Xun = X_unnormalize(Xn,Par_info)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% DDDDD  RRRR   EEEEE   AA   MM   MM PPPPP    AA   CCCCCC KK  KK   AA   GGGGGG EEEEE %%
%% DDDDDD RRRR   EEEEE  AAAA  MM   MM PPPPPP  AAAA  CCCCC  KK  KK  AAAA  GG     EEEEE %%
%% DD  DD RR RR  EE    AA  AA MMM MMM PP  PP AA  AA CC     KK KK  AA  AA GG     EE    %%
%% DD  DD RR RR  EEE   AA  AA MMMMMMM PP  PP AA  AA CC     KKKK   AA  AA GG  GG EEE   %%
%% DD  DD RRRRR  EEE   AAAAAA MMM MMM PPPPPP AAAAAA CC     KKKK   AAAAAA GG  GG EEE   %%
%% DD  DD RR RR  EE    AAAAAA MM   MM PPPPP  AAAAAA CC     KK KK  AAAAAA GG  GG EE    %%
%% DDDDDD RR  RR EEEEE AA  AA MM   MM PP     AA  AA CCCCC  KK  KK AA  AA GGGGGG EEEEE %%
%% DDDDD  RR  RR EEEEE AA  AA MM   MM PP     AA  AA CCCCCC KK  KK AA  AA GGGGGG EEEEE %%
%%                                                                                    %%
%% This function back transforms to original space the normalized parameter values    %%
%%                                                                                    %%
%% SYNOPSIS: Xun = X_unnormalize(Xn,Par_info)                                         %%
%% where                                                                              %%
%%  Xn          [input] Nxd matrix of normalized parameter values (N vectors)         %%
%%  Par_info    [input] Parameter structure: Ranges, initial/prior and bnd treatmnt   %%
%%  Xun         [outpt] Nxd matrix of parameter values in original space              %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Dec. 2016                                            %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

switch Par_info.norm
    case 0 % No normalization
        Xun = Xn;
    case 1
        % Latest MATLAB release
        % Xun = Par_info.minun + Xn * (Par_info.maxun - Par_info.minun);
        % Safer implementation
        Xun = bsxfun(@plus,bsxfun(@times,Xn,Par_info.maxun - ...
            Par_info.minun),Par_info.minun);
        % Safest implementation: should work with old MATLAB releases
        % [M,d] = size(Xun); Xun = nan(M,d);
        % for ii = 1:M
        %     Xun(ii,1:d) = Par_info.minun + Xn(ii,1:d) .* ...
        %         (Par_info.maxun - Par_info.minun);
        % end
end

end