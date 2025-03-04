function i = strmatch(str,strs,flag)
%STRMATCH Find possible matches for string.
%   I = STRMATCH(STR, STRARRAY) looks through the rows of the character
%   array or cell array of strings STRARRAY to find strings that begin
%   with the string contained in STR, and returns the matching row indices. 
%   Any trailing space characters in STR or STRARRAY are ignored when 
%   matching. STRMATCH is fastest when STRARRAY is a character array. 
%
%   I = STRMATCH(STR, STRARRAY, 'exact') compares STR with each row of
%   STRARRAY, looking for an exact match of the entire strings. Any 
%   trailing space characters in STR or STRARRAY are ignored when matching.
%
%   Examples
%     i = strmatch('max',strvcat('max','minimax','maximum'))
%   returns i = [1; 3] since rows 1 and 3 begin with 'max', and
%     i = strmatch('max',strvcat('max','minimax','maximum'),'exact')
%   returns i = 1, since only row 1 matches 'max' exactly.
%   
%   STRMATCH will be removed in a future release. Use STRNCMP instead.
%
%   See also STRFIND, STRVCAT, STRCMP, STRNCMP, REGEXP.

%   Mark W. Reichelt, 8-29-94
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.21.4.11 $  $Date: 2009/11/16 22:27:34 $

% The cell array implementation is in @cell/strmatch.m

if( nargin < 2 )
    error(nargchk(2,3,nargin,'struct'));
end

[m,n] = size(strs);
len = numel(str);

if (nargin==3)
    exactMatch = true;
    if ~ischar(flag)
        warning('MATLAB:strmatch:InvalidFlagType', ...
            ['The third argument to STRMATCH is not valid and will not be recognized in \n', ... 
            'a future version of MATLAB. Use string ''exact'' instead.']);
    elseif ~strcmpi(flag,'exact')
        warning('MATLAB:strmatch:InvalidFlag', ...
            ['Use flag ''exact'' in place of ''%s''. Flag ''%s'' \n', ...
            'is not valid and will not be recognized in a future version of MATLAB.'], flag, flag);
    end
else
    exactMatch = false;
end

% Special treatment for empty STR or STRS to avoid
% warnings and error below
if len==0
    str = reshape(str,1,len);
end 
if n==0
    strs = reshape(strs,max(m,1),n);
    [m,n] = size(strs);
end

if len > n
    i = [];
else
    if exactMatch && len < n % if 'exact' flag, pad str with blanks or nulls
        [strm,strn] = size(str);
        if strn ~= len
            error('MATLAB:strmatch:InvalidShape', ...
            'The first argument to STRMATCH should be a row vector when using flag ''exact''.');
        else
            % Use nulls if anything in the last column is a null.
            null = char(0); 
            space = ' ';
%            any(char(strs(:,end))==null)
            if ~isempty(strs) && any(char(strs(:,end))==null) 
                str = [str null(ones(1,n-len))];
            else
                str = [str space(ones(1,n-len))];
            end
            len = n;
        end
    end

    mask = true(m,1); 
    % walk from end of strs array and search for row starting with str.
    for outer = 1:m
        for inner = 1:len
            if (char(strs(outer,inner)) ~= str(inner))
                mask(outer) = false;
                break; % exit matching this row in strs with str.
            end   
        end
    end 
    i = find(mask);
end
