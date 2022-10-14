%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blkrltriang.m                                                 %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns for given matrices or systems A_1, ..., A_anz a 
% block right lower triangular matrix or system with blocks A_1, ...,A_anz.
% Therefore there must exist some N with anz = N * (N + 1) / 2. 
% As an example the output for anz = 6 is
%      /  0   0  A_1 \
% A = |   0  A_2 A_3  |
%      \ A_4 A_5 A_6 /
% Moreover, if A_i is not on the off diagonal of the resulting matrix, 
% A_i = 'z' results in a zero matrix of correct size and 
% A_i = 'e' results in an identity matrix of correct size.
%
% ----- Input ---------------------------------------------------------- 
%   varargin - List of objects such as matrices, LTI systems, ...
% ----- Output ---------------------------------------------------------
%          A - Object as described above
%
function [ A ] = blkrltriang( varargin )

anz = length(varargin);            % Number of inputs
n   = -0.5 + sqrt(0.25 + 2 * anz); % Block dimensions

% Sanity check
if n - floor(n) ~= 0
    error('blkrltriang:Wrong number of input arguments');
end

% Some initializations
k     = cell(n, 1);
m     = cell(n, 1);
sk    = cell(n + 1, 1);
sm    = cell(n + 1, 1);
sk{1} = 1;
issys = 0;

% Get the dimensions of the block matrices and sum them up by considering
% the off diagonal
temp  = 1;
for i = 1 : n
    [k{i}, ~] = size(varargin{temp});
    sk{i + 1} = sk{i} + k{i};
    temp      = temp + i;
end
temp  = temp - n;
sm{1} = 1;
for i = 1:n
   [~, m{n - i + 1}] = size(varargin{temp});
   sm{i + 1}         = sm{i} + m{n - i + 1};
   temp              = temp - n + i;
end

% Check whether an ss object is contained in the input.
for i = 1:anz
   if isa(varargin{i}, 'ss')
      issys = 1;
      break;
   end
end

% Construct the block right lower triangular matrix.
index = 1;
A     = zeros(sk{n + 1} - 1, sm{n + 1} - 1);
if issys
   A = ss(A); % Convert to ss if an ss object is involved
end
for i = 1 : n
    for j = 1 : n
       if i >= n - j + 1
          eci = sk{i + 1} - 1;
          eri = sm{j + 1} - 1;
          % Handle shortcuts
          if     isa(varargin{index}, 'char') && varargin{index} == 'z'
              A(sk{i}:eci, sm{j}:eri) = zeros(k{i}, m{n - j + 1});
          elseif isa(varargin{index}, 'char') && varargin{index} == 'e'
              A(sk{i}:eci, sm{j}:eri) = eye(k{i}, m{n - j + 1});
          else
              [t1, t2] = size(varargin{index});
              % Another sanity check
              if eci - sk{i} + 1 ~= t1 || eri - sm{j} + 1~= t2
                error('blkrltriang:check dimensions please')
              end
              A(sk{i}:eci, sm{j}:eri) = varargin{index};
          end
          index = index + 1;
       end
    end
end

end

