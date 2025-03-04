function str = num2words_rat(num,opts,varargin)
% Convert a numeric value to a string with an improper fraction in words. A NUM2WORDS example.
%
% (c) 2014-2020 Stephen Cobeldick
%
% Convert a numeric scalar decimal value into a string of words giving an
% improper fraction of the value. The numerator and denominator values are
% derived using the MATLAB function RAT, the tolerance may also be supplied.
% If the denominator is one, then the numerator only is returned.
%
%%% Syntax:
%  str = num2words_rat(num)
%  str = num2words_rat(num,opts)
%  str = num2words_rat(num,<name-value pairs>)
%
% File dependency: requires the function NUM2WORDS (FEX 47221).
%
%% Options %%
%
% The options may be supplied either
% 1) in a scalar structure, or
% 2) as a comma-separated list of name-value pairs.
%
% This function supports most of the same options as NUM2WORDS: see the
% function NUM2WORDS for the complete list. The changed options are:
%
% Field  | Permitted  |
% Name:  | Values:    | Description (example):
% =======|============|====================================================
% tol    | N          | Tolerance value N is passed directly to RAT.
% -------|------------|----------------------------------------------------
% type   |            | Is not supported by this function.
% -------|------------|----------------------------------------------------
%
%% Examples %%
%
% >> num2words_rat(0.25)
% ans = 'one fourth'
%
% >> num2words_rat(2.05, 'case','title', 'pos',true)
% ans = 'Positive Forty-One Twentieths'
%
% >> num2words_rat(2.3, 'tol',0.05)
% ans = 'seven thirds'
%
% % Mixed fraction:
% >> num = 13/5;
% >> int = num2words(fix(num));
% >> frc = num2words_rat(abs(rem(num,1)));
% >> str = sprintf('%s and %s', int, frc)
% str = 'two and three fifths'
%
%% Input and Output Arguments %%
%
%%% Inputs:
% num  = NumericScalar (float, int, or uint), the value to convert to words.
% opts = StructureScalar, with optional fields and values as per "Options" above.
% OR
% <name-value pairs> = a comma separated list of field names and associated values.
%
%%% Output:
% str = Char Vector, a string with the value of <num> in words of the chosen number
%       type and dialect, rounded to the requested order or significant figures.
%
% See also RAT NUM2WORDS WORDS2NUM NUM2SIP NUM2BIP NUM2ORDINAL INT2STR NUM2STR SPRINTF ARRAYFUN TTS

%% Input Wrangling %%
%
assert(isnumeric(num)&&isscalar(num),'SC:num2words_rat:NotScalarNumeric',...
	'First input <num> must be a numeric scalar.')
assert(isreal(num),'SC:num2words_rat:NotRealNumeric',...
	'First input <num> cannot be complex: %g%+gi',real(num),imag(num))
%
switch nargin
	case 1
		opts = struct();
	case 2
		assert(isstruct(opts)&&isscalar(opts),'SC:num2words_rat:NotScalarStruct',...
			'Second input <opts> can be a scalar structure.')
	otherwise
		opts = struct(opts,varargin{:});
		assert(isscalar(opts),'SC:num2words_rat:CellArrayValue',...
			'Invalid <name-value> pairs: cell array values are not permitted.')
end
%
fnm = fieldnames(opts);
%
assert(~any(strcmpi('type',fnm)),'SC:num2words_rat:TypeNotSupported',...
	'The option field name <type> is not supported')
%
%% Fraction and Numerator %%
%
idx = strcmpi('tol',fnm);
switch nnz(idx)
	case 0
		[N,D] = rat(num);
	case 1
		[N,D] = rat(num,opts.(fnm{idx}));
		opts = rmfield(opts,fnm{idx});
	otherwise
		error('SC:num2words_rat:DuplicateOption',...
			'Repeated field names:%s\b',sprintf(' <%s>,',fnm{idx}));
end
%
opts.type = 'decimal';
str = num2words(N,opts);
ino = +(1~=abs(N));
%
%% Denominator %%
%
switch D
	case 1
		return
	case 2
		vec = {'Half','Halve'};
		tmp = vec{1+ino};
		%	case 4
		%		tmp = 'Quarter';
	otherwise
		opts.pos = false;
		opts.type = 'ordinal';
		tmp = num2words(D,opts);
end
%
plu = 's';
tmp = [tmp,plu(1:ino)];
%
%% Character Case %%
%
idx = strcmpi('case',fnm);
%
if ~any(idx)
	str = sprintf('%s %s',str,lower(tmp));
else
	val = opts.(fnm{idx});
	switch lower(val);
		case 'upper'
			str = sprintf('%s %s',str,upper(tmp));
		case {'lower','sentence'}
			str = sprintf('%s %s',str,lower(tmp));
		case 'title'
			str = sprintf('%s %s',str,tmp);
		otherwise
			error('SC:num2words_rat:InvalidCase',...
				'Character <case> value ''%s'' is not supported',val)
	end
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%num2words_rat