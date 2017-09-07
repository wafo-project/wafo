function dcc = cc2dcc(param,cc,ddef)
% CC2DCC Discretize a cycle count.
%
% CALL:  dcc = cc2dcc(param,cc,ddef);
%
%        dcc   = a two column matrix with discrete classes.
%        param = the parameter matrix.
%        cc    = a two column matrix with the (continuous) cycle count.
%        ddef  = 1 causes peaks to be projected upwards and troughs 
%                  downwards to the closest discrete level (default).
%              = 0 causes peaks and troughs to be projected to 
%                  the closest discrete level.
%              =-1 causes peaks to be projected downwards and the 
%                  troughs upwards to the closest discrete level.
%
% Example:
%  x = load('sea.dat');
%  tp = dat2tp(x);
%  rfc = tp2rfc(tp);
%  param = [-2, 2, 41];
%  dcc = cc2dcc(param,rfc);
%  u = levels(param);
%  Frfc = dcc2cmat(dcc,param(3));
%  cmatplot(u,u,{Frfc}, 4);
%
%  close all;
%  
% See also  cc2cmat, dcc2cmat, dcc2cc

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 28-Jul-2000
%   Now correct upper and lower limits of discretization.
% Created by PJ (Par Johannesson) 01-Nov-1999
%   This is a new version of 'mkdisc' in WAT

% Check input arguments

ni = nargin;
no = nargout;
%error(nargchk(2,3,ni));
narginchk(2,3)
if ni<3 
  ddef=[];
end 

if isempty(ddef)
  ddef = 1;
end

dcc = cc;

% Make so that minima is in first column
II = find(cc(:,1)>cc(:,2));
if ~isempty(II)
  dcc(II,1) = cc(II,2);
  dcc(II,2) = cc(II,1);
end

% Make discretization

a=param(1); b=param(2); n=param(3);
delta = (b-a)/(n-1);        % Discretization step
dcc = (dcc-a)/delta + 1;

if ddef == 0
  dcc = min(max(round(dcc),1),n);
elseif ddef == +1
  dcc(:,1) = min(max(floor(dcc(:,1)),1),n-1);
  dcc(:,2) = min(max(ceil(dcc(:,2)),2),n);
elseif ddef == -1
  dcc(:,1) = min(max(ceil(dcc(:,1)),2),n);
  dcc(:,2) = min(max(floor(dcc(:,2)),1),n-1);
else
  error(['Undefined discretization definition, ddef = ' num2str(ddef)]);
end

% 
if ~isempty(II)
  dummy = dcc(II,1);
  dcc(II,1) = dcc(II,2);
  dcc(II,2) = dummy;
end
