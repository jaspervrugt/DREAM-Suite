function [nr,Km,Kstd,Kopt] = analyze_opt(eta,SSR,K)

% How many trials converge to same optimum?
[R,ii] = sort(SSR); d = size(eta,2);
% Sort eta's
eta_s = eta(ii,:); eta_opt = eta_s(1,1:d); nr = 1;
% Now check how many converged to same point
for i = 2:size(eta,1)
    df = sum(abs(eta_opt - eta_s(i,1:d)));
    if df < 1e-3
        nr = nr + 1;
    else
        break
    end
end
Km = mean(K); Kstd = std(K); Kopt = mean(K(ii(1:nr)));