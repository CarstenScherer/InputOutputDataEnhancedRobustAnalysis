%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : dgtrafo.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The dual multiplier corresponding to the standard DG scalings for a
% parametric uncertainty contained in the interval [a, b] and with oc
% occurences can be parametrized as
%          (*)^T * [0, M; M^T, 0] * [I, a*I; I, b*I] / (b-a)
% with M+M^T > 0. 
% This function returns the above outer factor in the case of multiple
% parametric uncertainties.
%
function [ F ] = dgtrafo( udata )
    % Some sanity checks
    arguments
        udata (:, 3) double
    end

    lu  = size(udata, 1);    % Number of uncertainty blocks
    
    F = [];
    for i = 1 : lu
        oc = udata(i, 3); % Number of occurrences/repetitions
        a  = udata(i, 1); % Start of parameter range
        b  = udata(i, 2); % End of parameter range

        % Outer factor for individual parametric uncertainty
        T = kron([1, a; 1, b], eye(oc));

        % Combination with previously obtained ones
        F = blkdiag(F, T);
    end
end