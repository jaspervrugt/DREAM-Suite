function Xn = X_normalize(Xun,Par_info)
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
%% This function transforms to normalized space the sampled parameter values          %%
%%                                                                                    %%
%% SYNOPSIS: Xn = X_normalize(Xun,Par_info)                                           %%
%% where                                                                              %%
%%  Xun         [input] Nxd matrix of parameter values in original space              %%
%%  Par_info    [input] Parameter structure: Ranges, initial/prior and bnd treatmnt   %%
%%  Xn          [outpt] Nxd matrix of normalized parameter values (N vectors)         %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Dec. 2016                                            %%
%% University of California Irvine                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

switch Par_info.norm
    case 0 % No normalization
        Xn = Xun;
    case 1
        % Latest MATLAB release
        % Xn = (Xun - Par_info.minun) ./ (Par_info.maxun - Par_info.minun);
        % Safer implementation
        Xn = bsxfun(@rdivide,Xun - Par_info.minun,Par_info.maxun - ...
            Par_info.minun);
        % Safest implementation: should work with old MATLAB releases
        % [M,d] = size(Xun); Xn = nan(M,d);
        % for ii = 1:M
        %     Xn(ii,1:d) = ( Xun(ii,1:d) - Par_info.minun ) ./ ...
        %         (Par_info.maxun - Par_info.minun);
        % end
end

end