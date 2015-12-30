function f = lh83pdf(t,h,mom,g)
%LH83PDF Longuet-Higgins (1983) approximation of the density (Tc,Ac) 
%         in a stationary Gaussian transform process X(t) where 
%         Y(t) = g(X(t)) (Y zero-mean Gaussian, X  non-Gaussian).
%         
%  CALL:  f   = lh83pdf(t,h,[m0,m1,m2],g);
%
%        f    = density of wave characteristics of half-wavelength
%               in a stationary Gaussian transformed process X(t),
%               where Y(t) = g(X(t)) (Y zero-mean Gaussian)
%       t,h   = vectors of periods and amplitudes, respectively.
%               default depending on the spectral moments
%    m0,m1,m2 = the 0'th,1'st and 2'nd moment of the spectral density
%               with angular  frequency.
%        g    = space transformation, Y(t)=g(X(t)), default: g is identity
%               transformation, i.e. X(t) = Y(t)  is Gaussian,
%               The transformation, g, can be estimated using lc2tr
%               or dat2tr or given apriori by ochi.
%       []    = default values are used.
%
% Example:
%  even = 0; nr =2;
%  mom  = spec2mom(jonswap,nr,[],even);
%  f    = lh83pdf([],[],mom);
%  pdfplot(f)
%
% See also  cav76pdf,  lc2tr, dat2tr

% References
% Longuet-Higgins,  M.S. (1983)
%"On the joint distribution wave periods and amplitudes in a 
% random wave field", Proc. R. Soc. A389, pp 24--258
%
% Longuet-Higgins,  M.S. (1975)
%"On the joint distribution wave periods and amplitudes of sea waves", 
% J. geophys. Res. 80, pp 2688--2694

% tested on: matlab 5.3 
% History:
% Revised pab 01.04.2001
% - Added example
% - Better automatic scaling for h,t
% revised by IR 18.06.2000, fixing transformation and transposing t and h to fit simpson req.
% revised by pab 28.09.1999
%   made more efficient calculation of f
% by Igor Rychlik


if nargin<3||isempty(mom)
  error('requires the moments')
elseif length(mom)~=3
  error('not enough moments')
else
  m0=mom(1);
  m1=mom(2);
  m2=mom(3);
end
 
if nargin<4||isempty(g)
   g=[(-5:0.02:5)' (-5:0.02:5)'];
   g(:,1)=g(:,1)*sqrt(m0);
end

L0=m0;
L1=m1/(2*pi);
L2=m2/(2*pi)^2;


if nargin<1||isempty(t)
    tt1=2*pi*sqrt(m0/m2);
    paramt=[0 1.7*tt1 51];
    t=levels(paramt);
end
    

if nargin<2||isempty(h) 
    px=gaus2dat([0. 0.;1 4.],g);
    px=abs(px(2,2)-px(1,2));
    paramh=[0 1.3*px 41];
    h=levels(paramh);
end


eps2 = sqrt((L2*L0)/(L1^2)-1);

if ~isreal(eps2)
  error('input moments are not correct')
end
const=4/sqrt(pi)/eps2/(1+1/sqrt(1+eps2^2));


a   = length(h); 
b   = length(t); 
der = ones(a,1);

h_lh = tranproc([h' der],g);
der  = abs(h_lh(:,2));
h_lh = h_lh(:,1);

% Normalization + transformation of t and h ???????
% Without any transformation

t_lh  = t/(L0/L1);
%h_lh = h_lh/sqrt(2*L0);
h_lh  = h_lh/sqrt(2);

t_lh  = 2*t_lh;


% Computation of the distribution
% if 0, %old call
%   for i=1:a,
%     for j=2:b,
%       f_th(i,j)=const*der(i)*(h_lh(i)/t_lh(j))^2*....
% 		exp(-h_lh(i)^2*(1+((1-1/t_lh(j))/eps2)^2))/((L0/L1)*sqrt(2)/2);
%     end
%   end  
% else %new call
  [T,H] = meshgrid(t_lh(2:b),h_lh);
  f_th=zeros(a,b);
  der=der(:);
  f_th(:,2:b)=const*der(:,ones(1,b-1)).*(H./T).^2.*....
      exp(-H.^2.*(1+((1-1./T)/eps2).^2))/((L0/L1)*sqrt(2)/2);
%  end

f      = createpdf(2);
f.f    = f_th;
f.x{1} = t(:);
f.x{2} = h(:);
f.labx = {'Tc','Ac'};
f.title = 'Joint density of (Tc,Ac) - Longuet-Higgins (1983)';
f.pl   = [10:20:90 95 99];
f.cl   = qlevels(f.f,f.pl);
%pdfplot(f)
