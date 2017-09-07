function [F_rfc,mu_rfc] = mctp2rfm(F,c_m)
% MCTP2RFM  Calculates the rainflow matrix for a MCTP.
%
% CALL:  [Frfc,mu_rfc] = mtp2rfc(F);
%        [Frfc,mu_rfc] = mtp2rfc(F,c_m);
%
% Frfc   = Rainflow Matrix (Rainflow Intensity)       [n,n]
% mu_rfc = Rainflow Counting Intensity                [n,n]
%
% F      = Cell array of min-Max and Max-min matrices {1,2}
% F{1,1} = min-Max matrix                             [n,n]
% F{1,2} = Max-min matrix                             [n,n]
% c_m    = Intensity of local minima
%
% If the matrix F{1,2} (F{1,2}=[]) is not given, then the process will
% be assumed to be time-reversible.
%
% Calculates the rainflow matrix for a Markov chain of
%   turning points.
%
% Example: 
%   param = [-1 1 32]; u = levels(param);
%   F = mktestmat(param,[-0.2 0.2],0.15,2);
%   Frfc = mctp2rfm({F []});
%   cmatplot(u,u,Frfc);
%
%   close all;
%
% See also  rfm2mctp, smctp2rfm, mctp2arfm, cmatplot, mc2rfm

% References
%  
%  P. Johannesson (1999):
%  Rainflow Analysis of Switching Markov Loads.
%  PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%  Lund Institute of Technology.
%  
%  P. Johannesson (1998):
%  Rainflow Cycles for Switching Processes with Markov Structure.
%  Probability in the Engineering and Informational Sciences, 
%  Vol. 12, No. 2, pp. 143-175.
  
% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1997
%   Copyright (c) 1997 by P�r Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

% Check input arguments

ni = nargin;
no = nargout;
%error(nargchk(1,2,ni));
narginchk(1,2)
if ni < 2
  c_m=[];
end

[F_rfc,mu_rfc] = smctp2rfm(1,F,c_m);

