function logL = logLcmat(Fobs,F,def)
% LOGLCMAT  log-Likelihood of cycle matrix.
%
% computes the log-Likelihood of an observed
%   cycle matrix Fobs which has the expected cycle matrix F.
%
% The log-Likelihood is
%   logL = C + D,
%   C = log(N!)-sum(log(N_ij!)),  D = sum(N_ij*log(g_ij))
%
% logL = logLcmat(Fobs,F,def)
%
% Fobs  = Observation of cycle matrix
% F     = Expected cycle matrix
% def   = 0: Don't compute constant part, logL=D, (default)
%         1: Compute constant part, logL=C+D

if nargin<3, def=0; end

F = flipud(F)';       % Convert to PJ-def
Fobs = flipud(Fobs)'; % Convert to PJ-def

n = length(F);
N = sum(sum(Fobs));

F = F / sum(sum(F));

FF = F(:);
FFobs = Fobs(:);
FI = find(FF>0 & FFobs>0);
FIobs = find(FFobs>0);

% Copmute constant part of log-Likelihood (part C)

if def == 1
  logL = sum(log(1:N));
  for k = 1:length(FIobs)
    logL = logL - sum(log( 1:FFobs(FIobs(k)) ));
  end
else
  logL = 0;
end

% Copmute log-Likelihood (part D)

logL = logL + sum(FFobs(FI).*log(FF(FI)));

