function S = ggamspec(wn,N,M)
%GGAMSPEC Generalized gamma spectrum in dimensionless form
%
% CALL S = ggamspec(wn,N,M)
%
%  S  = spectral values, same size as wn.
%  wn = normalized frequencies, w/wp.
%  N  = scalar defining the decay of the high frequency part. (default 5)
%  M  = scalar defining the spectral width around the peak. (default 4)
% 
% GGAMSPEC return the generalized gamma spectrum in non-
% dimensional form:
%      S = G0.*wn.^(-N).*exp(-B*wn.^(-M))  for wn > 0
%        = 0                              otherwise
%where
%  B  = N/M;
%  C  = (N-1)/M;
%  G0 = B^C*M/gamma(C), Normalizing factor related to Bretschneider form
%
% Note that N = 5, M = 4 corresponds to a normalized
% Bretschneider spectrum.
%
% Example
% wn = linspace(0,4);
% N = 6; M = 2;
% S = ggamspec(wn,N,M);
% plot(wn,S)
% assert(S(20:23),[0.705228792417578, 0.851577790974086,...
%                  0.974063777767012, 1.066965081263566], 1e-10)
%
%  close all
%
% See also mkbretschneider, mkjonswap, mktorsethaugen


% Reference
%  Torsethaugen, K. (2004)
%  "Simplified Double Peak Spectral Model for Ocean Waves"
%  In Proc. 14th ISOPE

% History
% revised pab april 2007
% -updated help header
% By pab jan 2007

%error(nargchk(1,3,nargin))
narginchk(1,3)
if nargin<2 || isempty(N)
  N = 5; % High frequency exponent
end

if nargin<3 || isempty(M)
  M = 4; % spectral width parameter.
end


S = zeros(size(wn));

%for w>0 % avoid division by zero
k = find(wn>0);
if any(k)
  B  = N/M;
  C  = (N-1)/M;

%     % A = Normalizing factor related to Bretschneider form
%     A    = B^C*M/gamma(C);
%     S(k) = A.*wn(k).^(-N).*exp(-B*wn(k).^(-M));

    logwn = log(wn(k));
    logA    = (C*log(B)+log(M)-gammaln(C));
    S(k)  = exp(logA-N.*logwn-B.*exp(-M.*logwn));

end
