%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : channel.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 02                                                            %
% Date    : 23.06.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a plant with inputs d_1, ..., d_m and outputs e_1, ..., e_n, this 
% function returns the corresponding plant with 
% inputs d_ic(1), ..., d_ic(end) and outputs e_oc(1), ..., e_oc(end).
%
% ----- Input ---------------------------------------------------------- 
%   p          - Struct
%   p.sys      - Given LTI system
%   p.inp      - Vector of input dimensions
%   p.out      - Vector of output dimensions
%   
%   ic         - Index of desired input channels
%   oc         - Index of desired output channels
% ----- Output ---------------------------------------------------------
%   p        - Same struct with modified elements
% 
function [p] = channel(p, varargin)

ic = varargin{1};

N = length(varargin); 
if N == 1
    oc = ic;
elseif N == 2
    oc = varargin{2};
else
    error('channel::Too many inputs');
end

li = mat2cell(1 : size(p.sys, 2), 1, p.inp);
li = li(ic);

lo = mat2cell(1 : size(p.sys, 1), 1, p.out);
lo = lo(oc);

p.sys = p.sys(cell2mat(lo), cell2mat(li));
p.inp = p.inp(ic);
p.out = p.out(oc);

end

