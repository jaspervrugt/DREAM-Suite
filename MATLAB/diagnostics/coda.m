function result = coda(draws,vnames,info,fid)

num_draws = length(draws);
if num_draws < 50
    error('coda: at least 50 draws are required');
end


if nargout == 1 % don't print, return a structure
    pflag = 1;
else
    pflag = 0; % print to the command window
end

if nargin == 4
    if ~isstruct(info)
        error('coda: must supply options as a structure');
    end
    nflag = 0;
    [vsize, ~] = size(vnames); % user may supply a blank vnames argument
    if vsize > 0
        nflag = 1; % we have variable names
    end
    fields = fieldnames(info);
    nf = length(fields);
    q = 0.025; r = 0.01; s = 0.95;
    p1 = 0.2; p2 = 0.5;
    for i=1:nf
        if strcmp(fields{i},'q')
            q = info.q;
        elseif strcmp(fields{i},'r')
            r = info.r;
        elseif strcmp(fields{i},'s')
            s = info.s;
        elseif strcmp(fields{i},'p1')
            p1 = info.p1;
        elseif strcmp(fields{i},'p2')
            p2 = info.p2;
        end
    end

elseif nargin == 3
    if ~isstruct(info)
        error('coda: must supply options as a structure');
    end
    nflag = 0;
    [vsize, ~] = size(vnames); 
    if vsize > 0
        nflag = 1; % we have variable names
    end
    fields = fieldnames(info);
    nf = length(fields);
    q = 0.025; r = 0.01; s = 0.95;
    p1 = 0.2; p2 = 0.5;
    fid = 1;
    for i=1:nf
        if strcmp(fields{i},'q')
            q = info.q;
        elseif strcmp(fields{i},'r')
            r = info.r;
        elseif strcmp(fields{i},'s')
            s = info.s;
        elseif strcmp(fields{i},'p1')
            p1 = info.p1;
        elseif strcmp(fields{i},'p2')
            p2 = info.p2;
        end
    end

elseif nargin == 2
    nflag = 1; % we have variable names
    q = 0.025; r = 0.01; s = 0.95; % set default values
    p1 = 0.2; p2 = 0.5;
    fid = 1;

elseif nargin == 1
    nflag = 0; % no variable names
    q = 0.025; r = 0.01; s = 0.95; % set default values
    p1 = 0.2; p2 = 0.5;
    fid = 1;

else
    error('Wrong # of arguments to coda');
end

result.q = q;
result.r = r;
result.s = s;

[ndraw,nvar] = size(draws);

if nflag == 0 % no variable names make some up
    Vname = [];
    for i=1 : nvar
        Vname{i} = str2mat(['variable ',num2str(i)]);
    end

elseif (nflag == 1) % the user supplied variable names
    Vname = [];
    [tst_n, ~] = size(vnames);
    if tst_n ~= nvar
        error('Wrong # of variable names in coda -- check vnames argument');
    end
    nmax = min(nsize,16); % truncate vnames to 16-characters
    for i=1:nvar
        Vname{i} = vnames(i,1:nmax);
    end
end % end of nflag issue

% =======> do SACF diagnostics
nlag = 25;
aout = zeros(25,nvar);
for i=1:nvar
    aout(:,i) = sacf(draws(:,i),nlag,1);
end

% pull out sacf's at 1,5,10,25
aprt = zeros(nvar,4);
aprt(:,1) = aout(1,:)';
aprt(:,2) = aout(5,:)';
aprt(:,3) = aout(10,:)';
aprt(:,4) = aout(25,:)';

% ========> do Raftery-Lewis diagnostics
rafout =  raftery(draws,q,r,s);
rout = zeros(nvar,5);
for i=1:nvar
    rout(i,1) = rafout(i).kthin;
    rout(i,2) = rafout(i).nburn;
    rout(i,3) = rafout(i).n;
    rout(i,4) = rafout(i).nmin;
    rout(i,5) = rafout(i).irl;
end

% =========> do Geweke diagnostics
geweke = momentg(draws);

% =========> split sample into 1st p1 percent and last p2 percent
%            and run Geweke chi-squared test
result(1).p1 = p1;
result(1).p2 = p2;
nobs1 = round(p1*ndraw);
nobs2 = round(p2*ndraw);
draws1 = draws(1:nobs1,:);
draws2 = trimr(draws,nobs2,0);

res1 = momentg(draws1);
res2 = momentg(draws2);

resapm = apm(res1,res2);

if pflag == 0 % print results to command window

    % =======> print SACF diagnostics
    fprintf(1,'MCMC CONVERGENCE diagnostics \n');
    fprintf(1,'Based on sample size = %10d \n',ndraw);
    fprintf(1,'Autocorrelations within each parameter chain \n');


    vstring = 'Variable';
    lstring1 = 'Lag 1';
    lstring2 = 'Lag 5';
    lstring3 = 'Lag 10';
    lstring4 = 'Lag 25';

    cnames = char(lstring1,lstring2,lstring3,lstring4);
    rnames = vstring;
    for i=1:nvar
        rnames = char(rnames,Vname{i});
    end
    in.fmt = '%12.3f';
    in.fid = fid;
    in.cnames = cnames;
    in.rnames = rnames;

    mprint(aprt,in);

    % print results with vnames
    fprintf(fid,'Raftery-Lewis Diagnostics for each parameter chain \n');
    fprintf(fid,'(q=%6.4f, r=%8.6f, s=%8.6f)\n',q,r,s);
    cstring1 = 'Thin';
    cstring2 = 'Burn';
    cstring3 = 'Total(N)';
    cstring4 = '(Nmin)';
    cstring5 = 'I-stat';

    cnames = char(cstring1,cstring2,cstring3,cstring4,cstring5);
    in2.fmt = char('%10d','%10d','%10d','%10d','%10.3f');
    in2.cnames = cnames;
    in2.fid = fid;
    in2.rnames = rnames;

    mprint(rout,in2);

    % =========> print Geweke diagnostics

    fprintf(fid,'Geweke Diagnostics for each parameter chain \n');
    cs1 = 'Mean';
    cs2 = 'std dev';
    cs3 = 'NSE iid';
    cs4 = 'RNE iid';

    cnames = char(cs1,cs2,cs3,cs4);
    in3.fmt = '%12.6f';
    gout = zeros(nvar,4);
    for i=1:nvar
        gout(i,:) = [geweke(i).pmean geweke(i).pstd geweke(i).nse geweke(i).rne];
    end
    in3.cnames = cnames;
    in3.rnames = rnames;
    in3.fid = fid;

    mprint(gout,in3);

    cs1 = 'NSE 4% ';
    cs2 = 'RNE 4% ';
    cs3 = 'NSE 8% ';
    cs4 = 'RNE 8% ';
    cs5 = 'NSE 15%';
    cs6 = 'RNE 15%';

    cnames = char(cs1,cs2,cs3,cs4,cs5,cs6);
    gout2 = zeros(nvar,6);
    for i=1:nvar
        gout2(i,:) = [geweke(i).nse1 geweke(i).rne1 geweke(i).nse2 geweke(i).rne2 ...
            geweke(i).nse3 geweke(i).rne3];
    end
    in4.cnames = cnames;
    in4.fid = fid;
    in4.fmt = '%12.6f';
    in4.rnames = rnames;
    mprint(gout2,in4);

    % =========> print Geweke chi-squared tests

    c = '%';
    % print results with vnames
    fprintf(1,'Geweke Chi-squared test for each parameter chain \n');
    fprintf(1,'First %2.0f%s versus Last %2.0f%s of the sample \n',100*p1,c,100*p2,c);
    clear in;
    in.cnames = char('Mean','N.S.E.','Chi-sq Prob');
    in.rnames = char('NSE estimate','i.i.d.','4% taper','8% taper','15% taper');
    in.fmt = '%12.6f';
    in.fid = fid;
    for i=1:nvar
        fprintf(1,'Variable %16s\n', strjust(Vname{i},'right'));
        gout3 = zeros(4,3);
        for k=1:4
            gout3(k,1) = resapm(i).pmean(k);
            gout3(k,2) = resapm(i).nse(k);
            gout3(k,3) = resapm(i).prob(k);
        end

        mprint(gout3,in);
    end

end % end of if pflag == 0

if pflag == 1 % return results structure


    result(1).nvar = nvar;
    result(1).meth = 'coda';
    result(1).ndraw = ndraw;

    for i=1:nvar
        result(i).kthin = rafout(i).kthin;
        result(i).nburn = rafout(i).nburn;
        result(i).n = rafout(i).n;
        result(i).nmin = rafout(i).nmin;
        result(i).irl = rafout(i).irl;
    end

    for i=1:nvar
        for j=1:4
            if j == 1
                result(i).auto1 = aprt(i,j);
            elseif j == 2
                result(i).auto5 = aprt(i,j);
            elseif j == 3
                result(i).auto10 = aprt(i,j);
            elseif j == 4
                result(i).auto50 = aprt(i,j);
            end
        end
    end

    for i=1:nvar
        result(i).nse1 = geweke(i).nse1;
        result(i).rne1 = geweke(i).rne1;
        result(i).nse2 = geweke(i).nse2;
        result(i).rne2 = geweke(i).rne2;
        result(i).nse3 = geweke(i).nse3;
        result(i).rne3 = geweke(i).rne3;
    end

    for i=1:nvar
        result(i).pmean = geweke(i).pmean;
        result(i).pstd  = geweke(i).pstd;
        result(i).nse   = geweke(i).nse;
        result(i).rne   = geweke(i).rne;
    end

end % end of if pflag == 1