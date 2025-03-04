function s = logsumexp(x, dim)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%LOGSUMEXP Returns log(sum(exp(x),dim)) while avoiding numerical          %
% underflow                                                               %
%                                                                         %
% dim = 1 for taking the mean along each column of x (default)            %
% dim = 2 for taking the mean along each row of x                         %
%                                                                         %
% © Written by Michael Chen (sth4nth@gmail.com)                           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin == 1
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end

end