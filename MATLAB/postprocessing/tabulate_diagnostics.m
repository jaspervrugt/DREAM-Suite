function tabulate_diagnostics(method,DREAMPar,chain,iloc,sndwch)
% Function tabulates some results of DREAM Package and returns posterior samples, Pars

%method = upper(method);

% Print header to screen
fid_d = fopen(strcat(method,'_diagnostics.txt'),'w');
% Now determine size for printing
evalstr_file = strcat(method,{' '},'diagnostics file');
ii_string = size(char(evalstr_file),2);
% Now diff
d_ii = 76 - ii_string - 2; d_half = d_ii/2; d1 = floor(d_half); d2 = d_ii-d1;
plot_str = strcat('=');
for i = 1:d1-1
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,{' '},evalstr_file,{' '});
for i = 1:d2-1
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,'\n');
fprintf(fid_d,char(plot_str));
switch sndwch 
    case 0
        fprintf(fid_d,'Sandwich correction is INACTIVE\n');
    case 1
        fprintf(fid_d,'Sandwich correction is ACTIVE\n');
end
fprintf(fid_d,'\n');
fprintf(fid_d,'TABLE 3: %s SINGLE CHAIN CONVERGENCE DIAGNOSTICS\n',method);

% Now calculate the convergence diagnostics for individual chains using the CODA toolbox
for j = 1:DREAMPar.N
    % First calculate diagnostics
    diagnostic_info{1} = coda(chain(floor(0.5*iloc):iloc,1:DREAMPar.d,j));
    % Now write to file DREAM.out
    diagnostic_info{1}(1).chain_number = j; prt_coda(diagnostic_info{1},[],fid_d);
end
d_ii = 76 - ii_string - 4; d_half = d_ii/2; d1 = floor(d_half); d2 = d_ii - d1;
plot_str = strcat('=');
for i = 1 : d1 - 2
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,{' '},strcat('End',{' '},evalstr_file,{' '}));
for i = 1 : d2 - 2
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,'\n');
fprintf(fid_d,char(plot_str));
% Now close the file again
fclose(fid_d);

% Now print to screen or not (not on unix/linux)
% if ispc || ismac, eval(char(strcat('edit',{' '},method,'_diagnostics.txt'))), end