function edges = autorule(x, minx, maxx, hardlimits)
xrange = maxx - minx;
if ~isempty(x) && (isinteger(x) || islogical(x) || isequal(round(x),x))...
        && xrange <= 50 && maxx <= flintmax(class(maxx))/2 ...
        && minx >= -flintmax(class(minx))/2
    edges = integerrule(x,minx,maxx,hardlimits,getmaxnumbins());
else
    edges = scottsrule(x,minx,maxx,hardlimits);
end
end

function edges = scottsrule(x, minx, maxx, hardlimits)
% Scott's normal reference rule
if ~isfloat(x)
    x = double(x);
end
binwidth = 3.5*std(x)/(numel(x)^(1/3));
if ~hardlimits
    edges = matlab.internal.math.binpicker(minx,maxx,[],binwidth);
else
    edges = matlab.internal.math.binpickerbl(min(x(:)),max(x(:)),minx,maxx,binwidth);
end
end