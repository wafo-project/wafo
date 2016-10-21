function [RFM,RFM0,res,RFMsid,RFM0sid,res_sid] = dtp2rfm_sid(x,y,n,ny,def)
%DTP2RFM_SID  Rainflow matrix from discrete turning points with side information.
%
% CALL:  [RFM,RFM0,res,RFM0sid,res_sid] = dtp2rfm(x,y,n,r,def)
%
% RFM     = Rainflow Matrix (residual included).    [n,n]
% RFM0    = Rainflow matrix (without resudual).     [n,n]
% res     = Residual.                               [2*n,1]
% RFMsid  = Rainflow Matrix with side information (residual included).    
%                                                   {r,r1}[n,n]
% RFMsid  = Rainflow Matrix with side information (without resudual).
%                                                   {r,r1}[n,n]
% res_sid = Residual for side information.          [2*n,1]
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
%   [RFM,RFM0,res,RFMsid,RFM0sid,res_sid] = dtp2rfm_sid(x,z,32,2,1);
%   figure(1),cmatplot(u,u,RFMsid,3)
%   RFM1 = RFMsid{1,1}+RFMsid{1,2}+RFMsid{2,1}+RFMsid{2,2};
%   figure(2),cmatplot(u,u,{RFM RFM1},3) % Shall be identical
%
% See also  dtp2arfm_sid, dtp2rfm, smctp2rfm

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

% Calculate asymmetric RFM and res

[RFM,RFM0,res,RFMsid,RFM0sid,res_sid] = dtp2arfm_sid(x,y,n,ny,def);

% Convert to symmetric rainflow

RFM = triu(RFM+RFM');
RFM0 = triu(RFM0+RFM0');
if def == 1
  for iy=1:ny
    for jy=1:ny
      RFMsid{iy,jy} = triu(RFMsid{iy,jy}+RFMsid{iy,jy}');
      RFM0sid{iy,jy} = triu(RFM0sid{iy,jy}+RFM0sid{iy,jy}');
    end
  end
end
if def == 2
  for iy=1:ny
    RFMsid{iy,1} = triu(RFMsid{iy,1}+RFMsid{iy,1}');
    RFM0sid{iy,1} = triu(RFM0sid{iy,1}+RFM0sid{iy,1}');
  end
end

