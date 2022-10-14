%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : test_multiple.m                                               %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Addpath for auxiliary functions.
addpath(genpath('../AuxiliaryFunctions'));

% Clean up.
clc
clear

%% Problem description

% Parameters involved in a satellite model essentially taken from
%  [G. F. Franklin, J. D. Powell, A. Emami-Naeini, Feedback Control of 
%   Dynamic Systems, 2010]
% We assume that k and b are unknown
J1 = 1;
J2 = 0.1;
k  = 0.091;          % Actual value of k that is only used for simulations
b  = 0.0036;         % Actual value of b that is only used for simulations
kI = [0.08, 0.12];   % Given interval containing k
bI = [0.0034, 0.02]; % Given interval containing b
% Corresponding uncertainties where we use the midpoint of the above
% intervals as nominal values.
ku    = ureal('ku', sum(kI)/2, 'Range', kI); 
bu    = ureal('bu', sum(bI)/2, 'Range', bI);
udata = [kI, 1; bI, 1];

% Two-degree-of-freedom design configuration with uncertainties pulled out.
% The interconnection signals are:
%   z1 = theta1 - theta2            w1 = ku * z1
%   z2 = d theta1 - d theta2        w2 = bu * z2
% For the states and remaining inputs and outputs see below.
A      = [0, 1, 0, 0; zeros(1, 4); 0, 0, 0, 1; zeros(1, 4)];
Bw     = [0, 0; [1, 1]/J2; 0, 0; [-1, -1]/J1];
Bnru   = [[0; 1; 0; 0], [0; 0; 0; 0], [0; 0; 0; 1/J1]];
Cz     = [-1, 0, 1, 0; 0, -1, 0, 1];
Dz     = zeros(2, 5);
Ce     = [1, 0, 0, 0; 0, 0, 0, 0];
Dewnru = [zeros(2, 3), [-1; 0], [0; 1]];
Cy     = [1, 0, 0, 0; 0, 0, 0, 0];
Dywnru = [zeros(2, 3), [0; 1], [0; 0]];

sys = ss(A, [Bw, Bnru], [Cz; Ce; Cy], [Dz; Dewnru; Dywnru], ...
         'StateName', {'theta2', 'd theta2', 'theta1', 'd theta1'}, ...
         'InputName', {'w1', 'w2', 'n', 'r', 'u'}, ...
         'OutputName', {'z1', 'z2', 'v-r', 'u', 'v', 'r'});

% Input and output partition
inp = [2, 1, 1, 1]; % [w, noise, reference, control] 
out = [2, 2, 2];    % [z, error, measurements]


% Discretize
Ts   = 0.05; % Sampling time
sysd = c2d(sys, Ts);

% Performance weights for the tracking error v-r, for the control u, for
% the noise n and the reference r, respectively.
s  = zpk('s');
we = c2d((0.5*s +0.433) / (s + 0.00433), Ts);
wu = 0.1;
wn = 0.4;
wr = 1;

% Weighted generalized plant (the uncertainty is still pulled out)
sysw = blkdiag(eye(out(1)), we, wu, eye(out(3))) * sysd * ...
       blkdiag(eye(inp(1)), wn, wr, eye(inp(4)));

% Weighted uncertain plant
syswu = lft(blkdiag(ku, bu), sysw);

% Design an Hinfty controller for the nominal weighted generalized plant
[K, ~, ~] = hinfsyn(syswu.NominalValue, out(3), inp(4));


%% Incorporating data

rng default; % For reproducibility

% First, we modify the weighted plant with uncertainties pulled out such
% that the measured signal y = [v; r] appears twice.
syswd = sysw([1:end, end-out(3)+1:end], :);
cl    = lft(syswd, K); % The last output signal of this is again y and this
                       % interconnection is now of the form as considered
                       % in the paper.

% Preparations for generating data from the actual system.
hor      = 30;              % Number of data points in each trajectory
timehor  = 0:Ts:(hor-1)*Ts; % Time-horizon of each trajectory

x0     = zeros(size(cl.a, 1), 1); % Known initial condition
cltrue = lft(blkdiag(k, b), cl);  % Actual weighted system

noisebnd = [0.1, 0.05, 0.01]; % Bound on the euklidian norm of the elements 
                              % of the noise sequence. Here, these elements 
                              % are scalar, which simplifies the generation 
                              % of the noise.
shor = [10, 15, 20, 25, 30];  % Some shorter horizons

ga = zeros(length(noisebnd), length(shor));

for l = 1 : length(noisebnd)
    % Generate data of the true system on full horizon
    n = noisebnd(l) * rand(hor, inp(2)); % Unknown noisy input
    r = 2*randn(hor, inp(3));            % Known input reference.
    d    = [n, r];                       % Generalized disturbance
    % We are only interested in the corresponding measured output.
    [y, ~, ~] = lsim(cltrue(out(2)+1:end, :), d, timehor, x0);

    
    for i = 1 : length(shor)
        yc{1}    = y(1:shor(i), :);
        rc{1}    = r(1:shor(i), :);
        toepcols = shor(i);
        tic
        ga(l, i) = ana_data(cl, udata, rc, yc, noisebnd(l), toepcols);
        toc
    end

end

