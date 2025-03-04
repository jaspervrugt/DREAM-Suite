function t_out = fix_time(t,n)

ii = find(isnan(t));
if isempty(ii)
    idx_ii = n;
else
    idx_ii = ii(1) - 1;
end
if idx_ii == 1
    t_out = t(1)*[1:n]';
else
    t_out = [t(1:idx_ii) ; t(idx_ii) + [1:(n-idx_ii)]'*(t(idx_ii)-t(idx_ii-1))];
end
