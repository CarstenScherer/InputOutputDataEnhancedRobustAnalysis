%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : toeplitzcols.m                                                %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given vector s = [s1; s2; ...; sh] where each sj is of length ls, 
% this function returns the first ncols columns of the Toeplitz matrix
%                       s1 0  0  ...  0
%                       s2 s1 0       :
%                       :   :    .    0
%                       sh sh-1  ... s1.
%
function [St] = toeplitzcols(s, ls, h, ncols)
    % Some sanity checks
    arguments
        s     (:, 1) double
        ls    (1, 1) {mustBeInteger}
        h     (1, 1) {mustBeInteger}
        ncols (1, 1) {mustBeInteger}
    end
    if h * ls ~= length(s)
        error('The input vector is of wrong dimension')
    end
    
    
    s  = [zeros(ls, 1); s];               % Prepend zeros
    s  = mat2cell(s, ls*ones(h+1, 1), 1); % Turn to cell
    T  = tril(toeplitz(1:h))+1; % Some indices for Toeplitz matrices
    T  = T(:, 1:ncols);         % We only take ncols columns
    St = cell2mat(s(T));        % Toeplitz matrix
end

