function [RFM,RFM0,res,RFMsid,RFM0sid,res_y] = dtp2arfm_sid(x,y,n,ny,def)
%DTP2ARFM_SID  Asymmetric RFM from discrete TP with side information.
%
% CALL:  [RFM,RFM0,res,RFMsid,RFM0sid,res_sid] = dtp2arfm_sid(x,y,n,r,def)
%
% RFM     = Rainflow Matrix (residual included).    [n,n]
% RFM0    = Rainflow matrix (without resudual).     [n,n]
% res     = Residual.                               [2*n,1]
% RFMsid  = Rainflow Matrix with side information (residual included).    
%                                                   {r,r1}[n,n]
% RFMsid  = Rainflow Matrix with side information (without resudual).
%                                                   {r,r1}[n,n]
% res_sid = Residual for side information.         [2*n,1]
%
% x     = Turning points (taking values 1,...,n).    [T,1]
% y     = Side information (taking values 1,...,r).  [T,1]
% n     = Number of levels.
% r     = Number of levels in side information..
% def   = Which type of side information
%          1: Mark min & max, r1=r
%          2: Mark when counted, r1=1
% RFM0  = Rainflow matrix (without resudual).     [n,n]
% res   = Residual.                               [2*n,1]
%
% Example: (Two processes as in Example 4.1 in PhD thesis)
%   P = [0.9 0.1; 0.05 0.95];
%   param = [-1 1 32]; u = levels(param);
%   F1 = mktestmat(param,[-0.4 -0.3],0.15,1);
%   F2 = mktestmat(param,[0.3 0.4],0.15,1);
%   [x,z] = smctpsim(P,{F1 F1'; F2 F2'},10000); % Two regime states
%   [RFM,RFM0,res,RFMsid,RFM0sid,res_sid] = dtp2arfm_sid(x,z,32,2,1);
%   figure(1),cmatplot(u,u,RFMsid,3)
%   RFM1 = RFMsid{1,1}+RFMsid{1,2}+RFMsid{2,1}+RFMsid{2,2};
%   figure(2),cmatplot(u,u,{RFM RFM1},3) % Shall be identical
%
% See also  dtp2arfm, smctp2arfm

% References:
%  
%  P. Johannesson (1999):
%  Rainflow Analysis of Switching Markov Loads.
%  PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%  Lund Institute of Technology.
  
% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  09-Apr-2001
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1998
% Copyright (c) 1997-1998 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.1, 22-Jan-1998

ni = nargin;
no = nargout;
error(nargchk(5,5,ni));

RFM0 = zeros(n);
res0 = [];
if def == 1
  RFM0sid = cell(ny,ny);
  for iy=1:ny
    for jy=1:ny
      RFM0sid{iy,jy} = RFM0;
    end
  end
  res0_y = [];
end
if def == 2
  RFM0sid = cell(ny,1);
  for iy=1:ny
    RFM0sid{iy,1} = RFM0;
  end
end


nres = length(res0);
res = zeros(2*n+1,1);
res(1:nres) = res0;
res_y = zeros(2*n+1,1);
res_y(1:nres) = res0_y;

% Calculate RFM and res

for k = 1:length(x)-1
  nres = nres+1;
  res(nres) = x(k);
  if def == 1
    res_y(nres) = y(k+1);
  end  
  cycleFound = 1;
  while cycleFound==1 & nres>=4
    A = sort([res(nres-1) res(nres-2)]);
    B = sort([res(nres) res(nres-3)]);
    if A(1) >= B(1) & A(2) <= B(2)
      i = res(nres-2);
      j = res(nres-1);
      RFM0(i,j) = RFM0(i,j) + 1;
      if def == 1
        s1 = res_y(nres-2);
        s2 = res_y(nres-1);
        RFM0sid{s1,s2}(i,j) = RFM0sid{s1,s2}(i,j) + 1;
        res_y(nres-2) = res_y(nres);
      end
      if def == 2
        s = y(k+1);
        RFM0sid{s,1}(i,j) = RFM0sid{s,1}(i,j) + 1;
      end

      res(nres-2) = res(nres);
      nres = nres-2;
    else
      cycleFound = 0;
    end
  end
end

% Residual

res = res(1:nres);
if def == 1
  res_y = res_y(1:nres);
end

% Calculate RFM = RFM0 + 'cycles in res'

RFM = RFM0;
if def == 1 | def == 2
  RFMsid = RFM0sid;
end

for k = 1:2:nres-1
  i = res(k);
  j = res(k+1);
  RFM(i,j) = RFM(i,j) + 1;
  if def == 1
    s1 = res_y(k);
    s2 = res_y(k+1);
    RFM0sid{s1,s2}(i,j) = RFM0sid{s1,s2}(i,j) + 1;
  end
  if def == 2
    s = y(end);
    RFM0sid{s,1}(i,j) = RFM0sid{s,1}(i,j) + 1;
  end
end


