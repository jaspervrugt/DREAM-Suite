function [pr_var,idcp] = extract_names(name)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Extract the variable names of univariate and multivariate priors      %%
%%                                                                       %%
%% SYNOPSIS: [pr_name,idp] = extract_names(name)                         %%
%%  where                                                                %%
%%   name      [input] Character string of prior handle                  %%
%%   pr_name   [outpt] Variable names of anomymous function handle prior %%
%%   idp       [outpt] Indices of closing parenthesis in character       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, June 2023                               %%
%% Los Alamos National Laboratory                                        %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

id = strfind(name,'('); id_s = id(1) + 1;   % Start of letters
idcp = strfind(name,')'); id_e = idcp(1)-1; % Index of closing parenthesis
name = name(id_s:id_e);                     % Variables separated by comma
id = strfind(name,','); n_id = numel(id);   % Number of variables
idp = nan(n_id,2);                          % 
for i = 1:n_id - 1
    idp(i,1:2) = [id(i) + 1 , id(i+1) - 1];
end
idp(n_id,1:2) = [id(n_id)+1 , numel(name)]; % Index of last name prior

pr_var = cell(1,size(idp,1));               % Initialize return argument
for i = 1:size(idp,1)
    pr_var{i} = name(idp(i,1):idp(i,2));    % Extract variable names handle
end

end