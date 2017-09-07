function amp_hist = cmat2amp(param,F)
%CMAT2AMP Calculates a histogram of amplitudes from a cycle matrix.
%
% CALL:  amp_hist = cmat2amp(param,F);
%
%   amp_hist = a two column matrix with amplitudes (defined by  param)
%              in the first column and frequencies in the second.
%   param    = the parameter matrix.
%   F        = the  nxn  frequency matrix for the cycle count.
%
% Example:
%   x = load('sea.dat');                   % Load data
%   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2
%   RFM = dtp2rfm(dtp,32);                 % Calculate rainflow matrix
%   amp_hist = cmat2amp([-2 2 32],RFM);    % Get amplitude histigram
%   bar(amp_hist(:,1),amp_hist(:,2));      % Plot histogram
%
%   close all;
%
% See also  cc2cmat 

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 03-Nov-1999

% Check input arguments

ni = nargin;
no = nargout;
%error(nargchk(2,2,ni));
narginchk(2,2)
n=param(3); % Number of discrete levels

amp_hist=zeros(n,2);


% First column: The values of the amplitudes
amp_hist(:,1) = levels([0 param(2)-param(1) n])'/2;

% Second  column: The number of amplitudes
for i=0:n-1
  amp_hist(i+1,2)=sum(diag(F,i));
end
