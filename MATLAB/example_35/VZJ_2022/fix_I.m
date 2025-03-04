function [idx_ii,ss] = fix_I(I,n,inf_rate,eta)

ii = find(isnan(I));
if isempty(ii)
    idx_ii = n;
else
    idx_ii = ii(1) - 1;
end
if abs(inf_rate(idx_ii) - eta(2)) < 1e-10
    ss = 1;
else
    ss = 0;
end