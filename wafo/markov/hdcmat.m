function hd = hdcmat(Fobs,F)
% HDCMAT  Hellinger distance of cycle matrix.
%
% comuptes the Hellinger distance between the observed
%   cycle matrix Fobs and the expected cycle matrix F.
%
% The Hellinger distance is defined as
%   H.D. = acos( sum( sqrt(N_ij/N*g_ij) ) )
%
% hd = hdcmat(Fobs,F)
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

hd = acos( sum( sqrt(FFobs(FI)/N.*FF(FI)) ) );

