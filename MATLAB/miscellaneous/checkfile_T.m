function T_new = checkfile_T(fname)
% Check content of restart budget file T.txt
% Written by Jasper A. Vrugt
% University of California Irvine
% DREAM Package

T_new = load(fname);
% Check the content of T_new
if isempty(T_new)       % --> empty
    error(['DREAM_PACKAGE ERROR: File ''T.txt'' is empty -->' ...
        ' Please store an integer in file ''T.txt'' ']);
end
if numel(T_new) > 1     % --> more than 1 variable listed
    error(['DREAM_PACKAGE ERROR: File ''T.txt'' lists more ' ...
        'than 1 variable --> Store only a single value ' ...
        '( = integer ) in file ''T.txt'' ']);
end
if ~isnumeric(T_new)    % --> not a numerical value listed
    error(['DREAM_PACKAGE ERROR: File ''T.txt'' does not store ' ...
        'a numerical value --> Store only a single value ' ...
        '( = integer ) in file ''T.txt'' ']);
end
if rem(T_new,1) > 0     % --> not an integer listed
    error(['DREAM_PACKAGE ERROR: File ''T.txt''stores a real' ...
        ' value with decimal fraction --> Store only a single' ...
        ' value ( = integer ) in file ''T.txt'' ']);
end
if (T_new <= 0)         % --> negative value
    error(['DREAM_PACKAGE ERROR: File ''T.txt'' stores ' ...
        'negative integers (or zero) -->  Store a single ' ...
        'positive integer in file ''T.txt'' ']);
    end

end
