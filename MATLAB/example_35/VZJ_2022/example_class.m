N = 100; clear x;
x = rand; fx = normpdf(x,0,1);
for i = 2:N
    z = x(i-1,1) + randn; 
    fz = normpdf(z,0,1);
    if fz > fx(i-1,1)
        x(i,1) = z; fx(i,1) = fz;
    else
        x(i,1) = x(i-1,1); fx(i,1) = fx(i-1,1);
    end
    % quick and dirty
%    F(i,1) = fx; X(i,1) = x; 
end