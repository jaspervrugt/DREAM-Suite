function [ Par_info ] = GPR_par_ranges ( func )

% Give the parameter ranges (minimum and maximum values)
Par_info.min(1) = 30 * 0.7696; Par_info.max(1) = 30 * 1.301;    % Corresponds to 1/0.05 to 1/0.17 ns/m in logarithmic units

% The x-axis is varying the fastest
x = 0 : 0.1 : 3;                        % Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
z = 0 : 0.1 : 3;                        % Boundaries of uniform z-grid

% Grid-cells in horizontal and vertical direction
func.dimhor = length(x) - 1; func.dimver = length(z) - 1;

% Scale DCT coefficients such that all models are possible
dum = zeros(func.dimver,func.dimhor);
count = 1;
for i = 1 : func.parz
    for j = 1 : func.parx
        if ( i > 1 || j > 1)
            count = count + 1;
            dum(i,j) = 1;
            dummy = idct2(dum);
            Par_info.min(count) = -1.7 / max(max(abs(dummy)));  % 0.2657
            Par_info.max(count) =  1.7 / max(max(abs(dummy)));  % 0.2657
            dum(i,j) = 0;
        end
    end
end

end