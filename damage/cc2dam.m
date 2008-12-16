function D=cc2dam(cc,beta,K)
% CC2DAM   Calculates the total Palmgren-Miner damage of a cycle count.
%
% CALL: D = cc2dam(cc,beta,K);
% 
%   D    = Damage.                                            [1xm]
%
%   cc   = Cycle count with minima in column 1 and            [nx2]
%          maxima in column 2.
%   beta = Beta-values, material parameter.                    [1xm]
%   K    = K-value, material parameter. (Optional, Default: 1) [1x1]
%
% The damage is calculated according to
%   D(i) = sum ( K * S^beta(i) ),  with  S = (max-min)/2
%
% Example:
%   x = load('sea.dat'); TP=dat2tp(x); RFC=tp2rfc(TP); 
%   bv = 3:8;
%   D = cc2dam(RFC,bv); plot(bv,D,'x-')
%
% See also  cmat2dam

% Tested on Matlab 6.0
%
% History:
% Revised by PJ  01-Nov-1999
% - updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.0'


% Check input and otput

ni = nargin;
no = nargout;
error(nargchk(2,3,ni));

if ni < 3
  K=[];
end

% Set default values

if isempty(K)
  K = 1;
end

% Calculate damage

amp = abs(cc(:,2)-cc(:,1))/2;

n=length(beta); D=zeros(1,n);
for i=1:n
  D(i)=K*sum(amp.^beta(i));
end

