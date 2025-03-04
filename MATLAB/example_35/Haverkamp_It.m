function y = Haverkamp_It(eta,plugin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Evaluates Haverkamp_I and Haverkamp_t and combines their infiltration times and    %%
%% corresponding cumulative infiltration values                                       %%
%%  SYNOPSIS: y = Haverkamp_It(eta,plugin)                                            %%
%% where                                                                              %%
%%  eta  [input]      4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]    %%
%%  plugin [input]    structure with Haverkamp input variables                        %%
%%  y    [outpt]      2nx1 vector of cum. infltr., I (cm), and time t (h)             %%
%%                                                                                    %%
%%  LITERATURE                                                                        %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               %%
%%      H. Vereecken (2023), The time validity of Philip's two-term infiltration      %%
%%      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       %%
%%      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration equation of    %%
%%      Parlange et al. (1982): Numerical solution, experimental design, and          %%
%%      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               %%
%%      https://doi.org/10.1002/vzj2.20167                                            %%
%%                                                                                    %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                                           %%
%%  University of California Irvine                                                   %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

I = Haverkamp_I(eta,plugin);        % First cumulative infiltration
t = Haverkamp_t(eta,plugin);        % Then time
y = [ t ; I ];                      % Return combined vector

end
