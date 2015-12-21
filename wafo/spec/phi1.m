function ph=phi1(w,h)
%PHI1 Factor transforming spectra to finite water depth spectra.
%
% CALL: tr = phi1(w,h)
%     
%    tr = vector of transformation factors
%     w = angular frequency
%     h = water depth
%
% Example: Transform a JONSWAP spectrum to a spectrum for waterdepth = 30 m
%   S = jonswap;
%   S1=S; S1.S=S1.S.*phi1(S1.w,30);
%
% 

% Reference
% Buows, E., Gunther, H., Rosenthal, W. and Vincent, C.L. (1985)
% 'Similarity of the wind wave spectrum in finite depth water: 1 spectral form.' 
%  J. Geophys. Res., Vol 90, No. C1, pp 975-986

% Tested on: Matlab 5.2
% History:
% by pab 16.02.2000
g=gravity;
if h==inf, % special case infinite water depth
  ph=ones(size(w));
  return
end
ph=zeros(size(w));

k1=w2k(w,0,inf);
dw1=2*w/g; % dw/dk|h=inf
k2=w2k(w,0,h);

dw2=ph;
ix=find(k1~=0);

dw2(ix)=dw1(ix)./(tanh(k2(ix)*h)+k2(ix)*h./cosh(k2(ix)*h).^2); % % dw/dk|h=h0
ph(ix)=(k1(ix)./k2(ix)).^3.*dw2(ix)./dw1(ix);