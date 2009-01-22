function [F_rfc,mu_rfc] = mc2rfm(Q,def)
% MC2RFM  Calculates the rainflow matrix/intensity for a Markov chain.
%
% CALL: [F_rfc,mu_rfc] = mc2rfc(Q,1);
%       [F_rfc,mu_rfc] = mc2rfc(Q,[2 h]);
%
% F_rfc  = Rainflow matrix / Rainflow intensity     [NxN]
% mu_rfc = Rainflow counting intensity              [NxN]
%
% Q      = Transition matrix for Markov chain       [nxn]
% def    = Definition 1: Markov chain  (default),    N=n
%                     2: Discretized Markov chain,   N=n+1
% h      = Discretization step (ONLY Def 2!)
%
% Calculates 
%   (1) the rainflow matrix for a Markov chain OR
%   (2) the rainflow intensity for a discretized Markov chain.
%
% Example: 
%   F = magic(5)
%   Q = mat2tmat(F)
%   Frfc = mc2rfm(Q)
%
% See also  smc2rfm, mctp2rfm, mc2stat, mc2reverse, cmatplot

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
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni<2, def = []; end
if isempty(def), def = 1; end


[F_rfc,mu_rfc] = smc2rfm(1,{Q},def);

