function [mM,Mm] = tp2mm(tp)
% TP2MM Calculates min2Max and Max2min cycles from turning points
%
% CALL:  [mM,Mm] = tp2mm(TP);
%
%   mM  = a two column matrix with the min2Max count.
%   Mm  = a two column matrix with the Max2min count.
%   TP  = a two column matrix with the sequence of turning points.
%
% Example:
%   x = load('sea.dat');
%   TP = dat2tp(x);
%   [mM,Mm] = tp2mm(TP);
%   ccplot(mM);
%
% See also  dat2tp, cc2cmat, ccplot

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 19-Oct-2000
%   Two versions existed (in 'onedim' and 'cycles')!
%   Removed version in 'onedim'
%   Now handles vectors
% Revised by PJ (Pär Johannesson) 01-Nov-1999
%   updated for WAFO
% Copied from WAT Ver. 1.2

[n m]= size(tp);
if n<m
  b=m;m=n;n=b; 
  tp=tp';
end

if n<2, 
  error('The vector must have more than 1 elements!')
end

switch m
  case {1, 2},  % dimension OK!
  otherwise, 
    error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ')
end

if tp(1,m)>tp(2,m)
  im = 2;
  iM = 1;
else
  im = 1;
  iM = 2;
end

% Delete first point if it is a maximum
%if tp(1,m)>tp(2,m)
%  tp = tp(2:n,:);
%  if tp(1,m)>tp(2,m)
%    error('tp  is not a sequence of turning points.')
%  end
%end

% Count min-max and max-min cycles
n=length(tp);
mM=[tp(im:2:n-1,m) tp(im+1:2:n,m)]; % min-max cycles
Mm=[tp(iM:2:n-1,m) tp(iM+1:2:n,m)]; % max-min cycles



