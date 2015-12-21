function [F_rfc,mu_rfc] = smc2rfm(P,Qc,def)
%SMC2RFM Calculates rainflow matrix/intensity for a switching Markov chain.
%
% [F_rfc,mu_rfc] = smc2rfc(P,Q,1);
% [F_rfc,mu_rfc] = smc2rfc(P,Q,[2 h]);
%
% F_rfc  = Rainflow matrix / Rainflow intensity     [NxN]
% mu_rfc = Rainflow counting intensity              [NxN]
%
% P      = Transition matrix for regime process     [rxr]
% Q      = Cell array of transition matrices        {r,1}
% Q{i}   = Transition matrix for Markov chain i     [nxn]
% def    = Definition 1: Markov chain (default),     N=n
%                     2: Discretized Markov chain,   N=n+1
% h      = Discretization step (ONLY Def 2!)
%
% Calculates 
%   (1) the rainflow matrix for a switching Markov chain OR
%   (2) the rainflow intensity for a discretized switching Markov chain.
%
% Example: 
%   P = [0.9, 0.1;0.05 0.95];
%   param = [-1 1 32]; u = levels(param);
%   [F1,F2] = mktestmat(param,[-0.3 0.3],0.15,1,-Inf);
%   Frfc = smc2rfm(P,{F1;F2});
%   cmatplot(u,u,Frfc)
%
% See also  mc2rfm, smctp2rfm, smc2stat, cmatplot

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
% Revised by PJ 19-May-2000
%   Changer disp(...) to warning(...) .
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997


% Check input arguments

ni = nargin;
%no = nargout;
error(nargchk(2,3,ni));

if ni<3, def = []; end
if isempty(def), def = 1; end

% Define 

Zstr = '123456789';

r = length(P);   % Number of regime states
n = length(Qc{1});  % Number of levels

% Check that the rowsums of P are equal to 1

sumP = sum(P');
if sum(sumP == 1) ~= length(P)
  warning('Rowsums of P not equal to 1. Renormalizing.');
  for i = 1:length(P)
    P(i,:) = P(i,:)/sumP(i);
  end
end

% Check that the rowsums of Qc{1},...,Qc{r} are equal to 1

for i = 1:r
  sumQi = sum(Qc{i}');
  if sum(sumQi == 1) ~= length(Qc{i})
    warning(['Rowsums of Q{' Zstr(i) '} not equal to 1. Renormalizing.']);
    for j = 1:length(Qc{i})
      Qc{i}(j,:) = Qc{i}(j,:)/sumQi(j);
    end
  end
end


% Make the transition matrix Q for the joint MC (X_k,Z_k)

Q = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  QQ = kron(Qc{z},P);
  Q(I+z,:) = QQ(I+z,:);
end


% Stationary distribution = ro of Q

ro = mc2stat(Q);

% Calculate rainflow matrix / rainflow intensity

if def(1) == 1  % Markov Chain

  N = n;
  mu_rfc = zeros(N,N);
  EYE = eye(n*r,n*r);
  Qcumsum = fliplr(cumsum(fliplr(Q)')');

  for i=2:n
    for j=i-1:n-1
      q = Qcumsum(:,r*j+1);

      A = Q(r*(i-1)+1:r*j,r*(i-1)+1:r*j); % i:j,   i:j
      C = Q(1:r*(i-1),r*(i-1)+1:r*j);     % 1:i-1, i:j
      Eye = EYE(1:r*(j-i+1),1:r*(j-i+1)); % eye (size of A)
      d = q(1:r*(i-1));                   % 1:i-1
      e = q(r*(i-1)+1:r*j);               % i:j
      Ro = ro(1:r*(i-1));                 % 1:i-1

      if j == i-1  % isempty(A)
        mu_rfc(i,j) = Ro*d;
      else
        mu_rfc(i,j) = Ro*(d+C*inv(Eye-A)*e);
      end
    end
  end

elseif def(1) == 2  % Discretized Markov Process

  N = n+1;
  mu_rfc = zeros(N,N);
  EYE = eye(n*r,n*r);
  Qcumsum = fliplr(cumsum(fliplr(Q)')');

  for i=2:N-1
    for j=i:N-1
      q = Qcumsum(:,r*(j-1)+1);

      A = Q(r*(i-1)+1:r*(j-1),r*(i-1)+1:r*(j-1)); % i:j-1, i:j-1
      C = Q(1:r*(i-1),r*(i-1)+1:r*(j-1));         % 1:i-1, i:j-1
      Eye = EYE(1:r*(j-i),1:r*(j-i));             % eye (size of A)
      d = q(1:r*(i-1));                           % 1:i-1
      e = q(r*(i-1)+1:r*(j-1));                   % i-1:j-1
      Ro = ro(1:r*(i-1));                         % 1:i-1

      if i == j  % isempty(A)
        mu_rfc(i,j) = Ro*d;
      else
        mu_rfc(i,j) = Ro*(d+C*inv(Eye-A)*e);
      end
    end
  end

end

% Fill the missing triangel in the mu_rfc matrix
if def(1) == 1

  F_rfc = zeros(n);
  for i = 1:n-1
    for j= i+1:n
      F_rfc(i,j) = mu_rfc(i+1,j-1)-mu_rfc(i,j-1)-mu_rfc(i+1,j)+mu_rfc(i,j);
    end
  end
  
  mu_rfc = cmat2nt(F_rfc);
  
% Old version
%  mu_rfc=flipud(mu_rfc'); % Convert to WAT matrix-format
%    for k = 2:N-1
%      for m = 1:N-k
%        i = k+m;
%        j = N+1+k-i;
%        mu_rfc(i,j) = mu_rfc(i,j-1)+mu_rfc(i-1,j)-mu_rfc(i-1,j-1);
%      end
%    end
%    
%  mu_rfc = flipud(mu_rfc)'; % Convert to PJ-def
%  
%  F_rfc = nt2cmat(mu_rfc);    % Calculate Rainflow matrix
%  
%%  F_rfc = fliplr(triu(fliplr(F_rfc),1)); % Zeros blow SW-NE diagonal

end

if def(1) == 2

  mu_rfc=flipud(mu_rfc'); % Convert to WAT matrix-format
  for k = 1:N-1
    for m = 1:N-k
      i = k+m;
      j = N+1+k-i;
      mu_rfc(i,j) = mu_rfc(i,j-1)+mu_rfc(i-1,j)-mu_rfc(i-1,j-1);
    end
  end

  % Calculate rainflow intensity by numerical differentiation

  h = def(2);
  F_rfc = zeros(N,N);

  for i=2:N-1
    for j=2:N-i
      F_rfc(i,j) = (mu_rfc(i-1,j-1)-mu_rfc(i-1,j+1) -(mu_rfc(i+1,j-1)-mu_rfc(i+1,j+1))) / (2*h)^2;
    end
  end
  
  F_rfc = flipud(F_rfc)';   % Convert to PJ-def
  mu_rfc = flipud(mu_rfc)'; % Convert to PJ-def

end



