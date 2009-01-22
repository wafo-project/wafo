function D = dspec2spreading(S)
%DSPEC2SPREADING Returns the directional spreading function D of S
%
%  CALL: D = dspec2spreading(S);
%
%  S = Directional spectrum structure with at least the fields:
%      .S     Spectrum values (size= [nt nf]).
%      .w OR .f  Frequency, length nf.
%      .theta Angular lags, length nt.
% D  = Directional spreading structure with the same fields as S.
%
% Let S(w,theta) = S(w)*D(w,theta), then D(w,theta) is the 
% directional spreading function of S, if int D(w,theta) dtheta = 1 for
% each w.
%
% See also  spreading


% Tested on: Matlab 6
% History
% by pab 13.10.2002
error(nargchk(1,1,nargin))

if (~isfield(S,'S'))
  error('S is missing from the spectrum structure!')
end

if (~isfield(S,'theta'))
  error('theta is missing from the spectrum structure!')
end

D   = S;
nt  = length(D.theta);
ind = find(D.S<0 | isnan(D.S));
if any(ind),
%  disp(['Negative directional distribution. Setting negative values to zero. min(DS) = '  num2str(min(DS(ind)))])
  D.S(ind) = 0;
end

%Normalize so that int D(theta,f) dtheta = 1 for each f 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Sf2      = simpson(D.theta,D.S);
Sf2       = trapz(D.theta,D.S);
k         = find(Sf2);
D.S(:,k)  = S.S(:,k)./Sf2(ones(nt,1),k);
%plot(Sf2)
%pause
return
