%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : ana.m                                                         %
%                                                                         %
% Author  : Tobias Holicki, University of Stuttgart                       %
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
% The computation of the upper bound relies on static multipliers and, more
% precisely, on static DG scalings. We employ here a dual version of the
% standard corresponding analysis criteria in order to compare these
% criteria with the other suggested ones.
% Finally, note that the use of static multipliers for constant parametric 
% uncertainties is known to be restrictive in general and one should use
% dynamic ones instead.
%
% ----- Input ---------------------------------------------------------- 
%   sys   - State-space model of the above known linear part.
%   udata - Information on the uncertainties. The j-th row of udata is
%            supposed to be [Intj(1), Intj(2), ocj].
%   opt   - Struct with fields
%       opt - Options used by the function mincx.
% ----- Output ---------------------------------------------------------
%   gao   - Obtained upper bound on the robust energy gain.
% 
function [ gao ] = ana(sys, udata, opt)
    % Some sanity checks
    arguments
        sys     {mustBeA(sys, "ss")}
        udata   (:, 3) double
        opt.opt (1, 5) double = [1e-3, 200, 1e7, 50, 1]
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
    
    %% Build outer factors
    
    [~, OY] = outerfactor(sys, inp, out, 'ana');
    
    % Transformation for nonnormalized DG-scalings
    F  = dgtrafo(udata);
    OY = blkdiag(eye(2*lx), F, eye(err+dis)) * OY;
    
    % Splitting for directly optimizing over ga
    OY1  = OY(1:end-dis, :); 
    OY2  = OY(end-dis+1:end, :);  
    
    %% Define variables
    
    setlmis([]);
    
    [ga, ~,  ~] = lmivar(1, [ 1, 1]); % Energy gain upper bound
    [ Y, ~, sY] = lmivar(1, [lx, 1]); % Lyapunov certificate
    
    % Static DG-scalings
    P = [];
    M = cell(lu, 1);
    for i = 1 : lu
        oc            = udata(i, 3); % Get occurrences/repetitions
        [M{i}, ~, sM] = lmivar(2, oc * [1, 1]);
        P             = blkdiag(P, [zeros(oc), sM; sM', zeros(oc)]);
    end
    
    % Inner matrices
    IY  = lmivar(3, blkdiag(sY, -sY, P, zeros(err)));
    IYr = blkdiag(zeros(2*lx + size(P, 2)), eye(err));
    
    %% Constraints
    
    % *Positivit LMI*
    k = newlmi;
    lmiterm([-k, 1, 1, Y], 1, 1);
    
    % *System LMI* (after a Schur complement)
    k = newlmi;
    lmiterm([-k, 1, 1, IY], OY1', OY1);
    lmiterm([-k, 1, 1,  0], OY1' * IYr * OY1);
    lmiterm([-k, 2, 1,  0], OY2);
    lmiterm([-k, 2, 2, ga], 1, eye(dis));
    
    % *Uncertainty LMIs*
    for i = 1 : lu
        k = newlmi;
        lmiterm([-k, 1, 1, M{i}], 1, 1, 's');
    end
    %% Solve the LMI problem
    
    lmis = getlmis;
    ndec = decnbr(lmis) - 1;
    
    [gao, ~] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);
    
    gao      = sqrt(gao);

end
