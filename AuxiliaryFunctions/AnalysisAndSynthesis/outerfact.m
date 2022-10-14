%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : outerfact.m                                                   %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 14.02.2018                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the standard outer factors appearing when dealing
% with dynamic multipliers
%
% ----- Input ---------------------------------------------------------- 
%   p         - struct
%   p.sys     - system
%   p.ins     - dimensions of input signals
%   p.out     - dimensions of output signals
%   p.ps      - filter psi for the primal LMI
%   p.ps2     - (optional) second filter for the primal LMI
%   p.ph      - (optional) filter phi for the dual LMI for synthesis
%   p.ph2     - (optional) second filter for the dual LMI
%
% ----- Output ---------------------------------------------------------
%   out       - struct
%   out.prim  - outer factor of primal LMI
%   out.dual  - outer factor of dual LMI (only if p.ph is used)
%
function [ out ] = outerfact( p )
%% Handle input data

% filters for the primal LMI
if ~isfield(p, 'ps2')
    p.ps2 = p.ps;
end

% filters for the dual LMI
if isfield(p, 'ph')
    if ~isfield(p, 'ph2')
        p.ph2 = p.ph;
    end
end


%% System and filter data
la   = size(p.sys.a, 1);  % number of system states
lio  = length(p.ins);     % number of input/output signals
lps1 = size(p.ps.a, 1);   % number of states in ps
lps2 = size(p.ps2.a, 1);  % number of states in ps2
lds1 = size(p.ps, 1);     % output dimension of ps
lds2 = size(p.ps2, 1);    % output dimension of ps2
if isfield(p, 'ph')
    lph1 = size(p.ph.a, 1);  % number of states in ph
    lph2 = size(p.ph2.a, 1); % number of states in ph2
    ldh1 = size(p.ph, 2);    % output dimension of ph^*
    ldh2 = size(p.ph2, 2);   % output dimension of ph2^*
end


% Realizations of filters
[Aps1, Bps1, Cps1, Dps1] = ssdata(p.ps);
[Aps2, Bps2, Cps2, Dps2] = ssdata(p.ps2);
if isfield(p, 'ph')
    [Aph1, Bph1, Cph1, Dph1] = ssdata(p.ph);
    [Aph2, Bph2, Cph2, Dph2] = ssdata(p.ph2);
end

% Realization of the system
[A, B, C, D] = sssdata(p.sys, p.out, p.ins);

%% Build outer factors
if isfield(p, 'ph') && lio == 2 % Synthesis with stability
    % Build annihilators
    U = blkdiag(eye(lps1+lps2), null([C{2}, D{2,1}]));
    V = blkdiag(eye(lph1+lph2), null([B{2}', D{1,2}']));
    
    % primal: realization of col(ps G, ps2)
    Ap = blkrutriang(Aps1, 'z', Bps1 * C{1}, Aps2, 'z', A);
    Bp = [Bps1 * D{1,1}; Bps2; B{1}];
    Cp = [Cps1, zeros(lds1, lps2), Dps1 * C{1}; ...
          zeros(lds2, lps1), Cps2, zeros(lds2, la)];
    Dp = [Dps1 * D{1,1}; Dps2];

    % dual: realization of col(ph^*, -ph2^* G^*)
    Ad = blkrutriang(-Aph1', 'z', 'z', -Aph2', -Cph2' * B{1}', -A');
    Bd = [Cph1'; -Cph2' * D{1,1}'; -C{1}'];
    Cd = [-Bph1', zeros(ldh1, la+lph2); ...
          zeros(ldh2, lps1), -Bph2', -Dph2' * B{1}'];
    Dd = [Dph1'; -Dph2' * D{1,1}'];
    
    % Outer Factors
    out.prim = [eye(size([Ap, Bp])); Ap, Bp; Cp, Dp] * U; % Primal LMI
    out.dual = [eye(size([Ad, Bd])); Ad, Bd; Cd, Dd] * V; % Dual LMI
    
elseif isfield(p, 'ph') && lio == 3 % Synthesis with performance
    % Build annihilators
    U = blkdiag(eye(lps1+lps2), null([C{3}, D{3,1}, D{3,2}]));
    V = blkdiag(eye(lph1+lph2), null([B{3}', D{1,3}', D{2,3}']));
    
    % primal: realization of [ps G, ps G12; ps2, 0; G21, G22; 0, I]
    Ap = blkrutriang(Aps1, 'z', Bps1 * C{1}, Aps2, 'z', A);
    Bp = [Bps1 * D{1,1}, Bps1 * D{1,2}; ...
          blklltriang(Bps2, B{1}, B{2})];
    Cp = blkrutriang(Cps1, 'z', Dps1*C{1}, Cps2, 'z', ...
                     [C{2}; zeros(p.ins(2), la)]); 
    Dp = [blklutriang(Dps1*D{1,1}, Dps1*D{1,2}, Dps2); ...
          blkrutriang(D{2,1}, D{2,2}, eye(p.ins(2)))];

    % dual: realization of [ph*, 0; -ph2*G*, -ph2*G21*; 0, I; -G12*, -G22*]
    Ad = blkrutriang(-Aph1', 'z', 'z', -Aph2', -Cph2' * B{1}', -A');
    Bd = [blklltriang(Cph1', -Cph2'*D{1,1}', -Cph2'*D{2,1}'); ...
          -C{1}', -C{2}'];
    Cd = blkrutriang(-Bph1', 'z', 'z', -Bph2', -Dph2' * B{1}', ...
                     [zeros(p.out(2), la); -B{2}']);
    Dd = [blklltriang(Dph1', -Dph2' * D{1,1}', -Dph2' * D{2, 1}'); ...
          blkrltriang(eye(p.out(2)), -D{1,2}', -D{2, 2}')];
    
    % Outer Factors
    out.prim = [eye(size([Ap, Bp])); Ap, Bp; Cp, Dp] * U;   % Primal LMI
    out.dual = [eye(size([Ad, Bd])); Ad, Bd; Cd, Dd] * V;   % Dual LMI
elseif ~isfield(p, 'ph') && lio == 1 % Stability 
    % primal: realization of col(ps G, ps2)
    Ap = blkrutriang(Aps1, 'z', Bps1 * C{1}, Aps2, 'z', A);
    Bp = [Bps1 * D{1,1}; Bps2; B{1}];
    Cp = [Cps1, zeros(lds1, lps2), Dps1 * C{1}; ...
          zeros(lds2, lps1), Cps2, zeros(lds2, la)];
    Dp = [Dps1 * D{1,1}; Dps2];
    
    % Outer Factors
    out.prim = [eye(size([Ap, Bp])); Ap, Bp; Cp, Dp];   % Primal LMI
    
elseif ~isfield(p, 'ph') && lio == 2 % Performance
    % primal: realization of [ps G, ps G12; ps2, 0; G21, G22; 0, I]
    Ap = blkrutriang(Aps1, 'z', Bps1 * C{1}, Aps2, 'z', A);
    Bp = [Bps1 * D{1,1}, Bps1 * D{1,2}; ...
          blklltriang(Bps2, B{1}, B{2})];
    Cp = blkrutriang(Cps1, 'z', Dps1*C{1}, Cps2, 'z', ...
                     [C{2}; zeros(p.ins(2), la)]); 
    Dp = [blklutriang(Dps1*D{1,1}, Dps1*D{1,2}, Dps2); ...
          blkrutriang(D{2,1}, D{2,2}, eye(p.ins(2)))];

    % Outer Factors
    out.prim = [eye(size([Ap, Bp])); Ap, Bp; Cp, Dp];   % Primal LMI
end


end

