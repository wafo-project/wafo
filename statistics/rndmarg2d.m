function [V,H,ind] = rndmarg2d(N,phat)
%RNDMARG2D Random points from a MARG2D distribution 
%
% CALL:  [R1,R2] = rndmarg2d(N,phat);
% 
%     R1,R2   = N random points in R^2.
%     N       = number of points generated
%     phat    = parameter structure array (see fitmarg2d)
%
%Example: Random points from a 2D Rayleigh distribution
%    x1=linspace(0,10)';
%    phat = createfdata('distribution',@pdfmarg2d,'params',[2 2 .5]);
%    phat.pdfoptions.distribution={'pdfray','pdfray'};
%    phat.pdfoptions.numpar =ones(1,2);
%    [y1,y2] = rndmarg2d(1000,phat);
%    f = pdfmarg2d(x1,x1,phat,'wdata',true,'mesh', true);
%    plot(f), hold on
%    plot(y1,y2,'.'), hold off
% 
%  See also  fitmarg2d , pdfmarg2d, cdfmarg2d


%   References:
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 p. 133-134.

%  tested on: matlab 5.2
% history
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vector to structure
%  Per A. Brodtkorb 28.01.99

if (nargin < 2), 
  error('Requires two input arguments.'); 
end

options = phat.pdfoptions;
PV = phat.params(1:options.numpar(1));
%PH = phat.params(options.numpar(1)+(1:options.numpar(2)));
%psi = phat.params{sum(options.numpar)+1};

try
  VDIST=getdistname(lower(options.distribution{1}));
  %HDIST=getdistname(lower(options.distribution{2}));
catch
  error('Not enough inputs')
end

% psi = options.psi; % interaction parameter
PV  = num2cell(PV(:).',1);
%PH  = num2cell(PH(:).',1);


rnd_v = str2func(['rnd',VDIST]);

rndsize = [N,1];
V = rnd_v(PV{:},rndsize);


% perform a direct inversion by newton
P=rand(rndsize);
[H ind]=invcmarg2d(V ,P,phat); % Inverse of the 2D  cdf given V . slow
end



