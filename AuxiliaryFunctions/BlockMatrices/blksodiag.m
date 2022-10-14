%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blksodiag.m                                                   %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 18.12.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns for given objects A_1, ... A_N the object
%
% A = diag([0, A_1'; A_1, 0], ... , [0, A_N'; A_N, 0])
%
% whenever this object is well-defined.
%
% ----- Input ---------------------------------------------------------- 
%   varargin - List of objects such as matrices, LTI systems, ...
%
% ----- Output ---------------------------------------------------------
%          A - An object as above
%
function [ A ] = blksodiag( varargin )

N = length(varargin); % Number of objects
A = cell(N, 1);       % Initialization

% Construct individual blocks
for i = 1 : N
    A{i} = blkodiag(varargin{i}', varargin{i});
end

% Block augmentation
A = blkdiag(A{:});

end

