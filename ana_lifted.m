%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : ana_lifted.m                                                  %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns (if possible) an upper bound on the robust energy
% gain of the uncertain feedback interconnection
%   x(t+1) =   Ax(t) +  Bw w(t) +  Bd d(t)
%     z(t) = Cz x(t) + Dzw w(t) + Dzd d(t)        w(t) = Del z(t)
%     e(t) = Ce x(t) + Dew w(t) + Ded d(t)
% where 
%           Del = blkdiag(del1 * I_oc1, ..., delN * I_ocN)
% and delj is a constant parameter which is known to be contained in some
% interval Intj.
%
% The computation of the upper bound relies on the discrete time lifting 
% procedure as explained, e.g., in [1] and on static multipliers. More
% precisely, we use static DG scalings and employ a dual version of the
% standard corresponding analysis criteria in order to compare these
% criteria with the other suggested ones.
%
% [1] T. Chen, B. A. Francis, Optimal Sampled-Data Control Systems, 1995.
%
% ----- Input ---------------------------------------------------------- 
%   sys   - State-space model of the above known linear part.
%   udata - Information on the uncertainties. The j-th row of udata is
%            supposed to be [Intj(1), Intj(2), ocj].
%   h     - Horizon for the lifting procedure.
%   opt   - Struct with fields
%       opt - Options used by the function mincx.
% ----- Output ---------------------------------------------------------
%   gao   - Obtained upper bound on the robust energy gain.
% 
function [ gao ] = ana_lifted( sys, udata, h, opt )
    % Some sanity checks
    arguments
        sys     {mustBeA(sys, "ss")}
        udata   (:, 3) double
        h       {mustBeInteger, mustBePositive} = 2
        opt.opt (1, 5) double   = [1e-3, 200, 1e7, 50, 1]
    end
    
    %% Abbreviations
    
    lx  = size(sys.a, 1);    % Dimension of system state
    lw  = sum(udata(:, 3));  % Dimension of uncertain signal
    err = size(sys, 1) - lw; % Dimension of error signal
    dis = size(sys, 2) - lw; % Dimension of generalized disturbance
    lu  = size(udata, 1);    % Number of uncertainty blocks

    % Partitioning of the relevant input and output signals
    inp = [udata(:, 3)', dis];
    out = [udata(:, 3)', err];
    
    % Generate corresponding lifted system
    [hsys, hinp, hout] = lifted_system(sys, inp, out, h);
    
    err = hout(end); % Size of lifted error signal
    dis = hinp(end); % Size of lifted generalized disturbance
    
    % Adjust uncertainty description.
    % We exploit here, that if the system is affected by an uncertainty
    % blkdiag(Del1, ..., DelN), then the lifted system is affected by the 
    % folloowing uncertainty involving a number of repetitions:
    %    blkdiag(kron(eye(h), Del1), ..., kron(eye(h), DelN)).
    udata(:, 3) = udata(:, 3) * h;
    
    %% Build outer factors
    
    [~, OY] = outerfactor(hsys, hinp, hout, 'ana');
    
    % Transformation for nonnormalized DG-scalings
    F  = dgtrafo(udata);
    OY = blkdiag(eye(2*lx), F, eye(err+dis)) * OY;
    
    % Splitting for directly optimizing over ga
    OY1  = OY(1:end-dis, :); 
    OY2  = OY(end-dis+1:end, :);  
    
    %% Define variables
    
    setlmis([]);
    
    [ga, ~,  ~] = lmivar(1, [1, 1]);  % Energy gain upper bound
    [ Y, n, sY] = lmivar(1, [lx, 1]); % Lyapunov certificate
    
    % Static DG-scalings for the lifted system. We actually use scalings
    % with a particular block band structure for better efficiency.
    P = [];
    M = cell(lu, 1);
    for i = 1 : lu
        oc            = udata(i, 3);
        % One could use [M{i}, ~, sM] = lmivar(2, oc*[1, 1]) but this only
        % yields slightly better upper bounds at the cost of a much larger
        % computational burden.
        [M{i}, n, sM] = blkbandlmivar(oc/h, h, n);
        P             = blkdiag(P, [zeros(oc), sM; sM', zeros(oc)]);
    end
    
    % Inner matrices
    IY  = lmivar(3, blkdiag(sY, -sY, P, zeros(err)));
    IYr = blkdiag(zeros(2*lx + size(P, 2)), eye(err));
    
    
    %% Constraints
    
    % *Positivit LMI*
    k = newlmi;
    lmiterm([-k, 1, 1, Y], 1, 1);
    
    % *Uncertainty LMIs*
    for i = 1 : lu
        k = newlmi;
        lmiterm([-k, 1, 1, M{i}], 1, 1, 's');
    end

    % *System LMI* (after a Schur complement)
    k = newlmi;
    lmiterm([-k, 1, 1, IY], OY1', OY1);
    lmiterm([-k, 1, 1,  0], OY1' * IYr * OY1);
    lmiterm([-k, 2, 1,  0], OY2);
    lmiterm([-k, 2, 2, ga], 1, eye(dis));
    

    
    %% Solve the LMI problem
    
    lmis = getlmis;
    ndec = decnbr(lmis) - 1;
    
    [gao, ~] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);
    gao      = sqrt(gao);

end