function [RFM,u,param] = dat2rfm(x,h,n)

%DAT2RFM  Calculates the rainflow matrix from a time signal.
%
% CALL:  RFM = dat2rfm(x)
%        [RFM,u,param] = dat2rfm(x,h,n)
%
% Input:
%   x   = Time signal.                          [N,1]/[N,2]
%   h   = Threshold range for rainflow filter. (default: 0)
%   n   = Number of discretization levels.     (default: 64) 
%         OR paramter matrix [a b n].
%
% Output:
%   RFM   = Rainflow matrix                       [N,N]
%   u     = Discrete levels.                      [n,1]
%   param = the parameter matrix  [a b n].
%
% Example:
%   x = load('sea.dat');
%   [RFM,u] = dat2rfm(x);    % Default parameters
%   subplot(1,2,1); cmatplot(u,u,RFM,3);
%   [RFM,u] = dat2rfm(x,0.5,[-2.5 2.5 50]);
%   subplot(1,2,2); cmatplot(u,u,RFM,3);
%
%   close all;
%
% See also  rfcfilter, dat2dtp, dtp2frm

% Copyright (c) 2003 by P�r Johannesson

% Tested  on Matlab  6.5
%
% History:
% Created by PJ (P�r Johannesson) 10-Apr-2003
% Updated by PJ 03-Jun-2003

%%%%
% Compile to MEX
%
% mcc -x ts2rfm

%%%%
% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,3,ni));

if ni<2, h=[]; end
if ni<3, n=[]; end

%%%%
% Default settings

if isempty(h), h=0; end
if isempty(n), n=64; end

if h<0, h=0; end

x = x(:,end);  % Data values in last column (skip time in 1st column)

%%%%
% Get Turning Points (TP)
% and 
% Rainflow filter signal

if h==0
    tp = rfcfilter(x,0,1);  % Get TP
else
    tp = rfcfilter(x,h);    % Get TP & RFC-filter
end

%%%%
% Discretization

if length(n)==3
    param = n;
    n=param(3);
else,
    u_min = min(tp);
    u_max = max(tp);

    param = [u_min u_max n];
end

u = levels(param);  % Discrete levels

%%%%
% Get discrete TP

[dtp,u] = dat2dtp(param,tp,h,0);

%%%%
% Calculate RFM

RFM = dtp2rfm(dtp,n,'CS');