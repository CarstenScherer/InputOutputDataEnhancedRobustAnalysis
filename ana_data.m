%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : ana_data.m                                                    %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns (if possible) an upper bound on the robust energy
% gain of the channel [n; r] -> e of the uncertain feedback interconnection
%   x(t+1) =   Ax(t) +  Bw w(t) +  Bn n(t) +  Br r(t)
%     z(t) = Cz x(t) + Dzw w(t) + Dzn n(t) + Dzr r(t)     w(t) = Del z(t)
%     e(t) = Ce x(t) + Dew w(t) + Den n(t) + Der r(t)
%     y(t) = Cy x(t) + Dyw w(t) + Dyn n(t) + Dyr r(t)
% where 
%           Del = blkdiag(del1 * I_oc1, ..., delN * I_ocN)
% and delj is a constant parameter which is known to be contained in some
% interval Intj.
%
% The computation of the upper bound relies on the discrete time lifting 
% procedure as explained, e.g., in [1] and on static multipliers. More
% precisely, we use a dual version of static DG scalings.
%
% The computation also relies on the availability of one or multiple noisy 
% input-output trajectories. For details see [2].
%
% [1] T. Chen, B. A. Francis, Optimal Sampled-Data Control Systems, 1995.
% [2] T. Holicki, C. W. Scherer, Input-Output-Data-Enhanced Robust Analysis
%     via Lifting, 2023.
%
% ----- Input ---------------------------------------------------------- 
% sys        - State-space model of the above known linear part.
% udata      - Information on the uncertainties. The j-th row of udata is
%              supposed to be [Intj(1), Intj(2), ocj].
% r          - Data of the reference signal, possibly of many trajectories.
% y          - Data of the measured signal, possibly of many trajectories.
% noisebound - Pointwise in time bound on the euclidian norm of the noise
%              signal.
% toepcols   - The data will be plugged into Toeplitz matrices. This number
%              determines how much columns of those matrices are utilized
%              in the end.
% opt        - Struct with fields
%              opt - Options used by the function mincx.
% ----- Output ---------------------------------------------------------
%   gao   - Obtained upper bound on the robust energy gain.
% 
function [ gao ] = ana_data( sys, udata, r, y, noisebound, toepcols, opt)
    % Some sanity checks
    arguments
        sys               {mustBeA(sys, "ss")}
        udata      (:, 3) double
        r                 cell
        y                 cell
        noisebound (1, 1) double
        toepcols   (1, 1) {mustBeInteger} = 1
        opt.opt    (1, 5) double          = [1e-3, 200, 1e11, 50, 1]
    end

    %% Abbreviations
    
    ntr = length(y);          % Number of measured trajectories
    h   = size(y{1}, 1);      % Horizon length of gathered data         
    ly  = size(y{1}, 2);      % Dimension of measured signal
    lr  = size(r{1}, 2);      % Dimension of reference signal
    lx  = size(sys.a, 1);     % Dimension of system state
    lw  = sum(udata(:, 3));   % Dimension of uncertain signal
    err = size(sys, 1)-lw-ly; % Dimension of error signal
    dis = size(sys, 2) - lw;  % Dimension of generalized disturbance
    ln  = dis - lr;           % Dimension of noise signal
    lu  = size(udata, 1);     % Number of uncertainty blocks


    % Partitioning of the relevant input and output signals
    inp = [udata(:, 3)', dis];
    out = [udata(:, 3)', err, ly];

    % Check wether the assumptions are satisfied
    [si, Mb, Md, N] = check_assumption(sys, inp, out, h);

    % Lifted system on reduced horizon
    [hsys, hinp, hout] = lifted_system(sys, inp, out, si);
    
    err = hout(end-1); % Size of lifted error signal
    dis = hinp(end);   % Size of lifted generalized disturbance
    
    % Adjust uncertainty description for the ones appearing in the lifted 
    % system.
    udatas = [udata(:, 1:2), udata(:, 3) * si];

    %% Build outer factors
    
    % The matrices corresponding to the measured output do not directly
    % appear in the outer factor, hence, they are not considered here.
    [~, OY] = outerfactor(hsys(1:sum(hout(1:end-1)), :), hinp, ...
                          hout(1:end-1), 'ana');

    % Transformation for nonnormalized DG-scalings
    F  = dgtrafo(udatas);
    OY = blkdiag(eye(2*lx), F, eye(err+dis)) * OY;

    % Splitting for directly optimizing over ga
    OY1  = OY(1:end-dis, :);
    OY2  = OY(end-dis+1:end, :);
    

    %% Build stuff for the data-based multiplier
    
    % Involved columns of Toeplitz matrices
    tY = [];
    tR = [];
    for i = 1 : ntr
        tY = [tY, toeplitzcols(reshape(y{i}', [], 1), ly, h, toepcols)];
        tR = [tR, toeplitzcols(reshape(r{i}', [], 1), lr, h, toepcols)];
    end
    
    % Lifted system on full horizon and relevant system matrices
    [hsys, hinp, hout] = lifted_system(sys, [udata(:, 3)', ln, lr], ...
                                       out, h);
    % Since we assume that trajectories are generated from zero initial
    % condition, we only need the corresponding D term.
    [~, ~, ~, D]       = sssdata(hsys, hinp, hout);

    % Build outer factor of the multiplier
    Dhzn = cell2mat(D(1:lu, lu+1));
    Dhzr = cell2mat(D(1:lu, lu+2));
    Dhyn = D{end, lu+1};
    Dhyr = D{end, lu+2};
    % Incorporate output and reference data
    OD = blkdiag([tY; tR], eye(ln*h))' * ...
         blkrltriang(eye(hout(end)), -[Dhzr'; Dhzn'], -[Dhyr'; Dhyn']);

    % This is the corresponding part in the main system LMI
    OD = OD * [[zeros(size(N, 2), lx); Mb'], blkdiag(N', Md')];

    % For the multiplier for Toeplitz matrix of noise, we use somewhat
    % unusual multipliers that exploit some of the Toeplitz structure, but
    % not all. They are closely related to static D-scalings
    ODs = cell(h, 1);
    for i = 1 : h
       ODs{i} = blkdiag(eye(toepcols*ntr), ...
                                    noisebound*sqrt(h-i+1)*eye(ln*h)) * OD;
    end

    %% Define variables
    
    setlmis([]);
    
    [ga, ~,  ~] = lmivar(1, [1, 1]);  % Energy gain upper bound
    [ Y, n, sY] = lmivar(1, [lx, 1]); % Lyapunov certificate

    % Static DG-scalings for the lifted system. We actually use scalings
    % with a particular block band structure for better efficiency.
    P = [];
    M = cell(lu, 1);
    for i = 1 : lu
        oc            = udatas(i, 3);
        % One could use [M{i}, ~, sM] = lmivar(2, oc*[1, 1]) but this only
        % yields slightly better upper bounds at the cost of a much larger
        % computational burden.
        [M{i}, n, sM] = blkbandlmivar(oc/si, si, n);
        P             = blkdiag(P, [zeros(oc), sM; sM', zeros(oc)]);
    end


    % Inner matrices
    IY  = lmivar(3, blkdiag(sY, -sY, P, zeros(err)));
    IYr = blkdiag(zeros(2*lx + size(P, 2)), eye(err));

    % Scalings for data part
    e    = @(k) diag(double(1:toepcols == k));
    entr = @(k) diag(double(1:ntr == k));
    q    = @(k) kron(diag(double(1:h >= k)), eye(ln));
    for j = 1 : ntr
        for i = 1 : toepcols
            [d(i, j), ~, sd] = lmivar(1, [1, 1]); 
            D{i, j} = lmivar(3, sd * blkdiag(kron(entr(j), e(i)), -q(i)));
        end
    end
    
    %% Constraints
    
    % *Positivit LMI*
    k = newlmi;
    lmiterm([-k, 1, 1, Y], 1, 1);
    
    % *Uncertainty LMIs* 
    % This corresponds to the employed DG-scalings
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
    % The above is the same as in the lifting approach. Below are the
    % modifications due to data.
    for j = 1 : ntr
        for i = 1 : toepcols
            lmiterm([-k, 1, 1, D{i, j}], ODs{i}', ODs{i});
        end
    end

    % *Data*
    for j = 1 : ntr
        for i = 1 : toepcols
            k = newlmi;
            lmiterm([-k, 1, 1, d(i, j)], 1, 1);
        end
    end

    %% Solve the LMI problem
    
    lmis = getlmis;
    ndec = decnbr(lmis) - 1;
    
    [gao, ~] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);
    gao      = sqrt(gao);

end

%% Auxilliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_assumptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [si, Mb, Md, N] = check_assumption(sys, inp, out, h)

    lx  = size(sys.a, 1);  % Number system states
    lu  = length(inp) - 1; % Number of uncertainty blocks
    le  = out(end-1);      % Dimension of error signal
    eps = 1e-9;            % Some small epsilon
    
    % Lifted system and its input and output dimensions
    [hsys, hinp, hout] = lifted_system(sys, inp, out, h);
    
    % We only need few of its describing matrices
    [~, B, ~, D] = sssdata(hsys, hinp, hout);
    Dhyw         = cell2mat(D(end, 1:lu));

    for si = h:-1:1
        Bsiw  = [];
        Dsiew = [];
        for j = 1 : lu
            lt    = (h-si)*inp(j); % An abbreviation

            % We generate lifted matrices on the shorter horizon si and 
            % extended those with zeros. These are then stacked
            % horizontally
            Bsiw  = [Bsiw, B{j}(:, lt+1:end), zeros(lx, lt)];       
            Dsiew = [Dsiew, D{lu+1, j}(1:si*le, lt+1:end), ...
                                                         zeros(si*le, lt)]; 

            % Construct the N matrices
            N{j}  = [eye(si*inp(j)), zeros(si*inp(j), lt)];
        end
        % Remaining essential matrices
        Mb = Bsiw * pinv(Dhyw);
        Md = Dsiew * pinv(Dhyw);
        
        % If the assumption on the kernels is satisfied for si, then
        % the following matrix differences should vanish.
        if norm(Mb * Dhyw - Bsiw) < eps && norm(Md * Dhyw - Dsiew) < eps
            break;
        end
        % The assumption is not satisfied if we did not find a suitable si.
        if si == 1
            error('The assumption on the kernels is not satisfied')
        end
    end
    N  = blkdiag(N{:});
end




