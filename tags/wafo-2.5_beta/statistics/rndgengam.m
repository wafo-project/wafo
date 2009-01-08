function R = rndgengam(varargin)
%RNDGENGAM Random matrices from a Generalized Gamma distribution.
%
% CALL:  R = rndgengam(a,b,c,sz);
%
%         R = matrix of random numbers
%    a,b,c  = parameters (see pdfgengam)
%      phat = Distribution parameter struct
%             as returned from FITGENGAM.  
%        sz = size(R)    (Default common size of a and b)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options). 
% Example:
%   R=rndgengam(1,2,4,1,100);
%   plotweib(R),shg
%
% See also pdfgengam, cdfgengam, invgengam, fitgengam, momgengam 

% Reference: 
% rejection sampling: there is no reference
% inversion method:
% Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% added ms 09.08.2000
% revised pab 23.10.2000
% %  - added default b and c
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R
%  - replaced inversion method with a modified version of 
%    rgamma from stixbox (Anders Holtsberg)
% The algorithm is a rejection method. The logarithm of the gamma 
% variable is simulated by dominating it with a double exponential.
% The proof is easy since the log density is convex!
% 
% Reference: There is no reference! Send me (Anders Holtsberg) an email
% if you can't  figure it out.

error(nargchk(1,inf,nargin))
Np = 3;
options = []; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});


if isempty(b), b=1;end
if isempty(c), c=1;end

if isempty(rndsize)
  [iscmn a b c] = iscomnsize(a,b,c);
else
  [iscmn a b c] = iscomnsize(a,b,c,zeros(rndsize{:}));
end
if ~iscmn
    error('a, b and c must be of common size or scalar.');
end

%
%R=invgengam(rand(size(a)),a,b,c); % slow
%return

R = zeros(size(a));

ok = ((a>0) & (b>0) & c>0);
k = find(ok);
if any(k),
  ak=a(k);
  y0 = log(ak)-1./sqrt(ak);
  c0  = ak - exp(y0);
  c1 =(ak.*y0 - exp(y0));
  c2 = abs((y0-log(ak)));
  
  accept=k;  omit=[];
  while ~isempty(accept)
    ak(omit)=[];  c0(omit) =[];
    c1(omit)=[];  c2(omit)=[];
    sz = size(ak);
    la = log(ak);
    y  = log(rand(sz)).*sign(rand(sz)-0.5)./c0 + la;

    f = ak.*y-exp(y) - c1;
    g = c0.*(c2 - abs(y-la));
  
    omit = find((log(rand(sz)) + g) <= f);
    if ~isempty(omit)
      R(accept(omit)) = exp(y(omit));
      accept(omit)=[];
    end
  end % while
  R(k) = R(k).^(1./c(k)).*b(k);
end
  
k1=find(~ok);
if any(k1);
  R(k1)=nan;
end
return


