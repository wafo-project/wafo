function kl = klcmat(Fobs,F)
% KLCMAT  Kullback-Leibler distance of cycle matrix.
%
% comuptes the Kullback-Leibler distance between the observed
%   cycle matrix Fobs and the expected cycle matrix F.
%
% The Kullback-Leibler distance is defined as
%   K.L. = sum( f_ij * log(N*f_ij/N_ij) )
%
% kl = klcmat(Fobs,F)
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
FI = find( (Fobs>0) & (F>0) );

% Compute Kullback-Leibler distance

kl = sum( FF(FI) .* log(N*F(FI)./Fobs(FI)) );

