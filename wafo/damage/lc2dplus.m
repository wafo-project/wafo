function D=lc2dplus(cross,beta,num_slice)
%LC2DPLUS Upper bound for total damage from level crossings.
%
%  Calculates the upper bound for the total damage/damage intensity using
%  the crossing spectrum/dowcrossing intensity.
%
%  CALL: D = lc2dplus(cross,beta,n);
%
%  where
%        D     = the upper bound for the total damage/damage intensity,
%        cross = a two column matrix with levels  u  and corresponding
%                downcrossing intensities,
%        beta  = a vector with beta-values,
%        n     = (optional input argument) the number of slice levels
%                (default = 500).

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

if nargin<3
  num_slice=500;
end

[cc,delta]=down2cc(cross,num_slice);

amplitudes=cc(:,1)-cc(:,2)/2;
cc_length=length(cc);

delta=max(cross(:,2))/num_slice;

deltas=[delta*ones(1,cc_length-1) .75*delta];
for i=1:length(beta)
  D(i)=deltas*amplitudes.^beta(i);
end
