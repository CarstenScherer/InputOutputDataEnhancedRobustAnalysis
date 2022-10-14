%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : zez.m                                                         %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 23.06.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For numbers k, l, m this functions returns the matrix
% [0_{k x l};
%  I_{l x l};
%  0_{m x l}].
%
% ----- Input ---------------------------------------------------------- 
%  varargin - The numbers k, l, m either directly or within a vector. 
%             If m is omitted it is set to zero.
% ----- Output ---------------------------------------------------------
%         M - The matrix as above.
%
function [M] = zez(varargin)

% Handle inputs
N = length(varargin);
if N == 1 && length(varargin{1}) == 2
    p = [varargin{1}, 0];
elseif N == 1 && length(varargin{1}) == 3
    p = varargin{1};
elseif N == 2
    p = [varargin{1}, varargin{2}, 0];
elseif N == 3
    p = [varargin{1}, varargin{2}, varargin{3}];
else
    error('zez::Check input please')
end

M = [zeros(p(1), p(2)); eye(p(2)); zeros(p(3), p(2))];

end

