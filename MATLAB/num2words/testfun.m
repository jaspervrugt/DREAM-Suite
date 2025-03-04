function chk = testfun(fun,swap,cmp)
% Test function for checking number<->word conversion functions.
%
% (c) 2011-2020 Stephen Cobeldick
%
% See also NUM2SIP_TEST NUM2BIP_TEST BIP2NUM_TEST SIP2NUM_TEST
% WORDS2NUM_TEST NUM2WORDS_TEST NUM2WORDSQ_TEST NUM2ORDINAL_TEST NUM2MYRIAD_TEST

tmp = {'%s %3d','<a href="matlab:opentoline(''%1$s'',%2$3d)">%1$s %2$3d</a>'};
fmt = tmp{1+feature('hotlinks')};
itr = 0;
cnt = 0;
wsp = @(s)regexprep(s,'\s',' ');
chk = @nestfun;
%
	function nestfun(varargin)
		% if swap:
		%    (out1, in2, in3,... fun, in1, out2, out3,...)
		% else:
		%    (in1, in2, in3,... fun, out1, out2, out3,...)
		%
		dbs = dbstack();
		%
		if ~nargin % post-processing
			fprintf('%s: %d of %d testcases failed.\n',dbs(2).file,cnt,itr)
			return
		end
		%
		itr = itr+1;
		%
		idx = find(cellfun(@(f)isequal(f,fun),varargin));
		assert(nnz(idx)==1,'Missing function handle.')
		%
		if swap
			varargin([1,idx+1]) = varargin([idx+1,1]);
		end
		%
		xfa = varargin(idx+1:end);
		ofa =  cell(size(xfa));
		boo = false(size(xfa));
		%
		[ofa{:}] = fun(varargin{1:idx-1});
		%
		for k = 1:numel(xfa)
			if ~isequal(class(ofa{k}),class(xfa{k}))
				boo(k) = true;
				otx = class(ofa{k});
				xtx = class(xfa{k});
			elseif ischar(xfa{k}) && ~cmp(wsp(ofa{k}),xfa{k})
				boo(k) = true;
				otx = sprintf('''%s''',ofa{k}.');
				xtx = sprintf('''%s''',xfa{k}.');
			elseif isnumeric(xfa{k}) && ~tfisequaln(ofa{k},xfa{k})
				boo(k) = true;
				otx = strrep(mat2str(ofa{k},23),' ',',');
				xtx = strrep(mat2str(xfa{k},23),' ',',');
			elseif iscellstr(xfa{k}) && ~isequal(wsp(ofa{k}),xfa{k})
				boo(k) = true;
				otx = sprintf(',''%s''',ofa{k}{:});
				xtx = sprintf(',''%s''',xfa{k}{:});
				otx = sprintf('{%s}',otx(2:end));
				xtx = sprintf('{%s}',xtx(2:end));
			end
			if boo(k)
				fprintf(fmt, dbs(2).file, dbs(2).line);
				fprintf(' (output argument %d)',k);
				fprintf('\noutput: %s\nexpect: %s\n', otx, xtx);
			end
		end
		cnt = cnt+any(boo);
	end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%testfun
function boo = tfisequaln(a,b)
% Compare equality of numeric matrices with 2*EPS tolerance.
af = isfinite(a(:));
bf = isfinite(b(:));
ai = isinf(a(:));
bi = isinf(b(:));
if isinteger(a) || isinteger(b)
	boo = isequal(a,b);
elseif isequal(size(a),size(b)) && ~any(xor(af,bf)|xor(ai,bi))
	boo = all(a(ai)==b(bi)) && all((abs(a(af)-b(bf))<=(eps(a(af))+eps(b(bf)))));
else
	boo = false;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tfisequaln