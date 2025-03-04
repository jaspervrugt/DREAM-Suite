function tabulate_output(method,DREAMPar,ML,MAP,MEAN,STD,CORR,sndwch)
% Function tabulates posterior moments and correlation coefficients of sampled chains

method = upper(method);
% Create legend/label string for different parameters
str_x0 = cell(DREAMPar.d,1);
for i = 1:DREAMPar.d
    str_x0(i,:) = cellstr(strcat('x_{',num2str(i),'}')); 
end

% Print header to screen
fid = fopen(strcat(method,'_output.txt'),'w');
% Now determine size for printing
evalstr_file = strcat(method,{' '},'output file'); ii_string = size(char(evalstr_file),2);
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
fprintf(fid,char(plot_str));
switch sndwch 
    case 0
        fprintf(fid,'Sandwich correction is INACTIVE\n');
    case 1
        fprintf(fid,'Sandwich correction is ACTIVE\n');
end
fprintf(fid,'\n');
fprintf(fid,'TABLE 1: %s POSTERIOR PARAMETER STATISTICS\n',upper(method));
fprintf(fid,'|-----------------------------------------------------|\n');
fprintf(fid,'| Parameter      ML        MAP      MEAN      STD     |\n');
fprintf(fid,'|-----------------------------------------------------|\n');

% Now print
for j = 1 : DREAMPar.d
    if (j < 10)
        fprintf(fid,'|   %s %+11.3f %+9.3f %+9.3f %+9.3f   |\n',char(str_x0(j)),ML(j),MAP(j),MEAN(j),STD(j));
    elseif (j >=10) && (j < 100)
        fprintf(fid,'|   %s %+10.3f %+9.3f %+9.3f %+9.3f   |\n',char(str_x0(j)),ML(j),MAP(j),MEAN(j),STD(j));
    else
        fprintf(fid,'|   %s %+9.3f %+9.3f %+9.3f %+9.3f   |\n',char(str_x0(j)),ML(j),MAP(j),MEAN(j),STD(j));
    end
end
fprintf(fid,'|-----------------------------------------------------|\n');
fprintf(fid,'ML: Maximum Likelihood values\n');
fprintf(fid,'MAP: Maximum A-Posteriori density values\n');
fprintf(fid,'MEAN: Posterior mean values\n');
fprintf(fid,'STD: Posterior standard deviation\n');
fprintf(fid,'\n');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Calculate correlation posterior samples
fprintf(fid,'TABLE 2: %s POSTERIOR CORRELATION COEFFICIENTS\n',method);
evalstr_top = strcat('fprintf(fid,''|-');
for j = 1 : (DREAMPar.d-1)*9 + 19
    evalstr_top = strcat(evalstr_top,'-');
end
evalstr_top = strcat(evalstr_top,'|\n'');');
eval(char(evalstr_top));

evalstr = strcat('fprintf(fid,'' %8.2f');
for k = 1 : DREAMPar.d - 1
    evalstr = strcat(evalstr,' %8.2f');
end
evalstr = strcat(evalstr,'   |\n''');
% Now add corr values
evalstr = strcat(evalstr,',corr_values(1:',num2str(DREAMPar.d),'));');
% Now print to screen

fprintf(fid,'|%9s',[   ]);
for i = 1 : DREAMPar.d-1
    fprintf(fid,'%9s',char(str_x0(i)));
end
fprintf(fid,'%9s  |',char(str_x0(DREAMPar.d)));

fprintf(fid,'\n');
for i = 1 : DREAMPar.d
    fprintf(fid,'|%8s',char(str_x0(i))); corr_values = CORR(i,1:DREAMPar.d);
    eval(char(evalstr));
end
eval(char(evalstr_top));
fprintf(fid,'\n');

d_ii = 76 - ii_string - 4; d_half = d_ii/2; d1 = floor(d_half); d2 = d_ii-d1;
plot_str = strcat('=');
for i = 1:d1-2
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,{' '},strcat('End',{' '},evalstr_file,{' '}));
for i = 1:d2-2
    plot_str = strcat(plot_str,'=');
end
plot_str = strcat(plot_str,'\n');
fprintf(fid,char(plot_str));
fclose(fid);

%if ( ispc || ismac ),  eval(char(strcat('edit',{' '},method,'_output.txt'))), end