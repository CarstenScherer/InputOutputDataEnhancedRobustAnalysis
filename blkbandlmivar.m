%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blkbandlmivar.m                                               %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a nonsymmtric LMI variable with block band structure.
% Each block is square of dimension blksize and we have numdiagonalblocks
% on the diagonal of this variable.
%
function [M, n, sM] = blkbandlmivar(blksize, numdiagonalblocks, ...
                                    numcurrentvars)
    % Some sanity checks
    arguments
        blksize           (1, 1) {mustBeInteger, mustBePositive}
        numdiagonalblocks (1, 1) {mustBeInteger, mustBePositive}
        numcurrentvars    (1, 1) {mustBeInteger, mustBePositive}
    end

    % Abbreviation
    h = numdiagonalblocks;

    % Build a matrix indicating the desired band structure
    I = diag(ones(h, 1)) + diag(ones(h-1, 1), 1) + diag(ones(h-1, 1), -1);
    I = kron(I, ones(blksize));

    % Now generate a matrix for use within LMILab
    numnewvars = length(find(I == 1)); % Number variables that we will add
    newvars    = numcurrentvars + (1:numnewvars);
    I(I == 1)  = newvars;

    % Generate new structured variable
    [M, n, sM] = lmivar(3, I);
end

