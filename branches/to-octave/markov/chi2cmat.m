function chi2 = chi2cmat(Fobs,F)
% CHI2CMAT  Chi-square distance of cycle matrix.
%
% Computes the chi-square distance between the observed
% cycle matrix Fobs and the expected cycle matrix F.
%
% The chi-square distance is defined as
%   chi2 = sum( (N_ij - N*g_ij)^2 ./ (N*g_ij) );
%
% chi2 = chi2cmat(Fobs,F)
%
% Fobs  = Observation of cycle matrix
% F     = Expected cycle matrix

F = flipud(F)';       % Convert to PJ-def
Fobs = flipud(Fobs)'; % Convert to PJ-def

n = length(F);
N = sum(sum(Fobs));

F = F / sum(sum(F));

FF = F(:);
FFobs = Fobs(:);
FI = find(F>0);

% Compute chi-square quantity

chi2 = sum( (FFobs(FI) - N*FF(FI)).^2 ./ (N*FF(FI)) );

