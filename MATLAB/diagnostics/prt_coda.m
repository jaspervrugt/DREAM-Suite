function prt_coda(results,vnames,fid)
% PURPOSE: Prints convergence diagnostics from DREAM samples
%---------------------------------------------------
% USAGE: prt_coda(results,vnames,fid)
% Where: results = a structure returned by coda, raftery, apm, momentg
%        fid     = file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%---------------------------------------------------
%  NOTES:   e.g. vnames = ['beta    ',
%                          'sigma   '];  NOTE: fixed width
%           e.g. fid = fopen('gibbs.out','wr');
%  use prt_coda(results,[],fid) to print to a file with no vnames
% --------------------------------------------------
%  RETURNS: nothing, just prints the diagnostic results
% --------------------------------------------------
% SEE ALSO: prt
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

if ~isstruct(results)
    error('prt_coda requires structure argument');
elseif nargin == 1
    nflag = 0; fid = 1;
elseif nargin == 2
    fid = 1; nflag = 1;
elseif nargin == 3
    nflag = 0;
    [vsize , ~] = size(vnames); % user may supply a blank argument
    if vsize > 0
        nflag = 1;
    end
else
    error('Wrong # of arguments to prt_coda');
end

nvar = results(1).nvar;

if nflag == 0 % no variable names make some up
    for i=1:nvar
%        Vname{i} = char(['variable ',num2str(i)]);
        Vname{i} = char(num2str(i));
    end

elseif (nflag == 1) % the user supplied variable names
    [tst_n, nsize] = size(vnames);
    if tst_n ~= nvar
        error('Wrong # of variable names in coda -- check vnames argument');
    end
    nmax = min(nsize,16); % truncate vnames to 16-characters
    for i=1:nvar
        Vname{i} = vnames(i,1:nmax);
    end
end % end of nflag issue

switch results(1).meth

    case {'raftery'} % <=================== raftery diagnostics

        % ========> Raftery-Lewis diagnostics
        rout = zeros(nvar,5);
        for i=1:nvar
            rout(i,1) = results(i).kthin;
            rout(i,2) = results(i).nburn;
            rout(i,3) = results(i).n;
            rout(i,4) = results(i).nmin;
            rout(i,5) = results(i).irl;
        end

        % print results with vnames
        fprintf(fid,'2. RAFTERY-LEWIS DIAGNOSTIC (q = %4.2f, r = %4.2f, s = %4.2f)\n',results(1).q,results(1).r,results(1).s);
        vstring  = 'Variable';
        cstring1 = 'Thin';
        cstring2 = 'Burn';
        cstring3 = 'Total(N)';
        cstring4 = '(Nmin)';
        cstring5 = 'I-stat';

        in.cnames = char(cstring1,cstring2,cstring3,cstring4,cstring5);
        in.fmt = char('%8d','%8d','%8d','%8d','%8.3f');
        in.fid = fid;
        rnames = vstring;
        for i=1:nvar
            rnames = char(rnames,Vname{i});
        end
        in.rnames = rnames;
        mprint(rout,in);

    case {'momentg'} % <=================== Geweke NSE, RNE diagnostics

        fprintf(fid,'3. GEWEKE DIAGNOSTIC\n');
        vs  = 'Variable    ';
        cs1 = 'Mean';
        cs2 = 'std dev';
        cs3 = 'NSE iid';
        cs4 = 'RNE iid';

        in0.cnames = char(cs1,cs2,cs3,cs4);
        in0.fmt = '%8.4f';
        gout = zeros(nvar,4);
        for i=1:nvar
            gout(i,:) = [results(i).pmean results(i).pstd results(i).nse results(i).rne];
        end

        rnames = vs;
        for i=1:nvar
            rnames = char(rnames,Vname{i});
        end

        in0.rnames = rnames;
        in0.fid = fid;

        mprint(gout,in0);

        cs1 = 'NSE 4% ';
        cs2 = 'RNE 4% ';
        cs3 = 'NSE 8% ';
        cs4 = 'RNE 8% ';
        cs5 = 'NSE 15%';
        cs6 = 'RNE 15%';

        in1.cnames = char(cs1,cs2,cs3,cs4,cs5,cs6);
        gout2 = zeros(nvar,6);
        for i=1:nvar
            gout2(i,:) = [results(i).nse1 results(i).rne1 results(i).nse2 results(i).rne2 ...
                results(i).nse3 results(i).rne3];
        end
        in1.fid = fid;
        in1.fmt = '%8.4f';
        in1.rnames = rnames;
        mprint(gout2,in1);

    case {'apm'} % <=================== Geweke chi-sqr diagnostics

        p1 = results(1).p1;
        p2 = results(1).p2;
        c = '%';
        fprintf(fid,'Geweke Chi-squared test for each parameter \n');
        fprintf(fid,'based on %d draws \n',results(1).ndraw);
        fprintf(fid,'First %2.0f%s versus Last %2.0f%s of the sample \n',100*p1,c,100*p2,c);

        in3.cnames = char('Mean','N.S.E.','Chi-sq Prob');
        in3.rnames = char('NSE estimate','i.i.d.','4% taper','8% taper','15% taper');
        in3.fmt = '%12.6f';
        for i=1:nvar
            fprintf(1,'Variable %16s\n', strjust(Vname{i},'right'));
            gout3 = zeros(4,3);
            for k=1:4
                gout3(k,1) = results(i).pmean(k);
                gout3(k,2) = results(i).nse(k);
                gout3(k,3) = results(i).prob(k);
            end
            in3.fid = fid;
            mprint(gout3,in3);
        end

    case{'coda'}
        nvar = results(1).nvar;
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,' CHAIN %d USING %d SAMPLES\n',results(1).chain_number,results(1).ndraw);
        fprintf(fid,'\n');
        fprintf(fid,'A. AUTOCORRELATION \n');
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
        in.fmt = '%10.3f';
        in.rnames = rnames;
        in.cnames = cnames;

        aprt = zeros(nvar,4);
        for i=1:nvar
            aprt(i,1) = results(i).auto1;
            aprt(i,2) = results(i).auto5;
            aprt(i,3) = results(i).auto10;
            aprt(i,4) = results(i).auto50;
        end
        in.fid = fid;
        mprint(aprt,in);

        fprintf(fid,'B. RAFTERY-LEWIS DIAGNOSTIC (q = %4.2f, r = %4.2f, s = %4.2f)\n',results(1).q,results(1).r,results(1).s);
        
        cstring1 = ' Thin ';
        cstring2 = ' Burn ';
        cstring3 = ' N     ';
        cstring4 = ' Nmin  ';
        cstring5 = 'I-stat';

        in2.cnames = char(cstring1,cstring2,cstring3,cstring4,cstring5);
        rout = zeros(nvar,5);
        for i=1:nvar
            rout(i,1) = results(i).kthin;
            rout(i,2) = results(i).nburn;
            rout(i,3) = results(i).n;
            rout(i,4) = results(i).nmin;
            rout(i,5) = results(i).irl;
        end
        fmt = char('%10d','%10d','%10d','%10d','%10.3f');
        %fmt = char('%10f'); % char('%10f','%10f','%10f','%10f','%10.3f');
        in2.fid = fid;
        in2.rnames = rnames;
        in2.fmt = fmt;
        mprint(rout,in2);

        % =========> print Geweke diagnostics

        fprintf(fid,'C. GEWEKE DIAGNOSTIC\n');
        vs  = 'Variable';
        cs1 = 'Mean';
        cs2 = 'std dev';
        cs3 = 'NSE iid';
        cs4 = 'RNE iid';

        in.cnames = char(cs1,cs2,cs3,cs4);
        in.fmt = '%10.3f';
        gout = zeros(nvar,4);
        for i=1:nvar
            gout(i,:) = [results(i).pmean results(i).pstd results(i).nse results(i).rne];
        end
        in.fid = fid;
        in.rnames = rnames;
        mprint(gout,in);

        cs1 = 'NSE 4% ';
        cs2 = 'RNE 4% ';
        cs3 = 'NSE 8% ';
        cs4 = 'RNE 8% ';
        cs5 = 'NSE 15%';
        cs6 = 'RNE 15%';

        in2.cnames = char(cs1,cs2,cs3,cs4,cs5,cs6);
        gout2 = zeros(nvar,6);
        for i=1:nvar
            gout2(i,:) = [results(i).nse1 results(i).rne1 results(i).nse2 results(i).rne2 ...
                results(i).nse3 results(i).rne3];
        end
        in2.fid = fid;
        in2.fmt = '%10.3f';
        in2.rnames = rnames;
        mprint(gout2,in2);

    otherwise
        error('results structure not known by prt_coda function');

end
