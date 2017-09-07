function [dtp,res]=rfm2dtp(F,res,N)
%RFM2DTP  Reconstructs a sequence of turning points from a rainflow matrix.
%
% CALL:  tp=rfm2dtp(RFM,res)
%        tp=rfm2dtp(RFM)
%        tp=rfm2dtp(RFM,res,N)
%
% Input:
%   RFM   = Rainflow matrix                        [n,n]
%   res   = Residual.                              [nres,1]
%   N     = Generate approximately N points
%
% Output:
%   tp    = Turning points.                        [N,1]
%
% Generates a sequence of turning points from a rainflow matrix and its
% residual.  If the residual is given, the cycles in the residual should
% not be included in the rainflow matrix.   The rainflow count of the
% output will be exactly the input rainflow matrix.
% If a residual is not given, then a stationary residual is generated.  
% With the third argument you can set the length of the output signal.
%
% Example:
%   x=load('sea.dat');
%   param = [-2 2 64]; n=param(3);
%   dtp0 = dat2dtp(param,x(:,2));
%   [RFM,RFM0,res0] = dtp2rfm(dtp0,n);
%   dtp = rfm2dtp(RFM0,res0);
%   plot(1:length(dtp0),dtp0,'b',1:length(dtp),dtp,'r')
%
% See also  dat2dtp, dtp2rfm

% Copyright (c) 2004 by Pär Johannesson

% Tested  on Matlab  6.5
%
% History:
% Created by PJ (Pär Johannesson) 16-Feb-2004

%%%%
% Check input arguments

ni = nargin;
%no = nargout;
%error(nargchk(1,3,ni));
narginchk(1,3)
if ni<2, res=[]; end
if ni<3, N=[]; end

f = fliplr(F)'; % Convert to FAT-def
n=length(F);
res = n-res+1;  % Convert to FAT-def

% Call function 'rfc2load' originally from FAT (Fatigue Analysis Toolbox)
% FAT is a predecessor of WAFO
[dtp,res] = rfc2load_fat(f,res,N);

dtp=n-dtp+1;    % Convert to WAFO-def
res=n-res+1;    % Convert to WAFO-def