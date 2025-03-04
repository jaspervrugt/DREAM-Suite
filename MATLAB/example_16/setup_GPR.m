function func = setup_GPR(func)
% Set-up forward kernel and create synthetic data

% The x-axis is varying the fastest
x = 0 : 0.1 : 3;                        % Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
z = 0 : 0.1 : 3;                        % Boundaries of uniform z-grid

sourcex = 0.01;                         % x-position of source
sourcez = 0.05 : 0.1 : 3.;              % z-positions of sources
receiverx = 2.99;                       % x-position of receivers
receiverz = 0.05 : .1 : 3;              % z-positions of receivers
nsource = length(sourcez); nreceiver = length(receiverz);

% Calculate acquisition geometry (multiple-offset gather)
for j = 1 : nsource
    for i = 1 : nreceiver
        data( ( j - 1 ) * nreceiver + i , :) = [sourcex sourcez(j) receiverx receiverz(i)];
    end
end

% Calculate forward modeling kernel, Courtesy Dr. James Irving, UNIL
func.J = tomokernel_straight(data,x,z); % Distance of ray-segment in each cell for each ray

% Grid-cells in horizontal and vertical direction
func.dimhor = length(x) - 1; func.dimver = length(z) - 1;

func.por = 0.36;                        % Porosity field
nz = 30;                                % Original discretization of water saturation model
nx = 40;
wcon = load('Sw.dat');                  % Load water saturation model: Courtesy of Dr. Michael Kowalsky (LBNL), Kowalsky et al. (2005; WRR)
count = 0;
for k = 1 : nz
    for i = 6 : nx - 5                  % Make model 3 by 3 m
        wcont( k , i - 5) = wcon( nz - k + 1 , i);
    end
end
wcont = wcont * func.por;               % Transform saturation into water content
wcont = dct2(wcont);                    % Transform water content model into DCT space

wtrunc = zeros(func.dimver,func.dimhor);
for j = 1 : func.parz
    for i = 1 : func.parx
        wtrunc(j,i) = wcont(j,i);       % Truncate at same order as inverse parameterization
    end
end
wtrunc = idct2(wtrunc);                 % Do inverse DCT
wtrunc = wtrunc';

% Translate water content model into slowness using the Refractive Index/CRIM model
func.pw = 81;                           % Permittivity of water
func.pa = 1;                            % Permittivity of air
func.ps = 5;                            % Permittivity of mineral grains
slowtrue = wtrunc(:) * sqrt(func.pw) + (func.por-wtrunc(:)) * sqrt(func.pa) + (1 - func.por) * sqrt(func.ps); % sqrt of effective permittivity
slowtrue = slowtrue / 0.3;              % True slowness field

% Simulate data for true slowness model
% Add Gaussian uncorrelated noise with a standard deviation of 1 ns.
func.datasim = func.J * slowtrue + func.error * randn(nsource * nreceiver,1); % Simulate data

end
