function f = cav76pdf(t,h,mom,g)
%CAV76PDF Cavanie et al. (1976) approximation of the density  (Tc,Ac)
%         in a stationary Gaussian transform process X(t) where 
%         Y(t) = g(X(t)) (Y zero-mean Gaussian, X  non-Gaussian).
%          
% CALL:  f = cav76pdf(t,h,[m0,m2,m4],g);
%
%        f    = density of wave characteristics of half-wavelength
%               in a stationary Gaussian transformed process X(t),
%               where Y(t) = g(X(t)) (Y zero-mean Gaussian)
%       t,h   = vectors of periods and amplitudes, respectively.
%               default depending on the spectral moments
% m0,m2,m4    = the 0'th, 2'nd and 4'th moment of the spectral density
%               with angular frequency.
%        g    = space transformation, Y(t)=g(X(t)), default: g is identity
%               transformation, i.e. X(t) = Y(t)  is Gaussian,
%               The transformation, g, can be estimated using lc2tr
%               or dat2tr or given a priori by ochi.
%       []    = default values are used.
%
% Example: 
%  mom = spec2mom(jonswap,4); 
%  f = cav76pdf([],[],mom);
%
% See also  lh83pdf, lc2tr, dat2tr

% References:
% Cavanie, A., Arhan, M. and Ezraty, R. (1976)
% "A statistical relationship between individual heights and periods of
%  storm waves".
% In Proceedings Conference on Behaviour of Offshore Structures,
% Trondheim, pp. 354--360
% Norwegian Institute of Technology, Trondheim, Norway
%
% Lindgren, G. and Rychlik, I. (1982)
% Wave Characteristics Distributions for Gaussian Waves --
% Wave-lenght, Amplitude and Steepness, Ocean Engng vol 9, pp. 411-432.

% tested on: matlab 5.3 NB! note
% History:
% revised pab 04.11.2000
% - fixed xlabels i.e. f.labx={'Tc','Ac'}
% revised by IR 4 X 2000. fixed transform and normalisation
% using Lindgren & Rychlik (1982) paper.
% At the end of the function there is a text with derivation of the density.
%
% revised by jr 21.02.2000
% - Introduced cell array for f.x for use with pdfplot 
% by pab 28.09.1999

if nargin<3||isempty(mom)
  error('requires the moments')
elseif length(mom)<3
  error('not enough moments')
else
  m0=mom(1);
  m2=mom(2);
  m4=mom(3);
end
 
if nargin<4||isempty(g)
  g=[(-5:0.02:5)' (-5:0.02:5)'];
  g(:,1)=sqrt(m0)*g(:,1);
end

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

eps4 = 1-m2^2/(m0*m4);
alfa = m2/sqrt(m0*m4);
if ~isreal(sqrt(eps4)) 
  error('input moments are not correct')
end



a   = length(h); 
b   = length(t); 
der = ones(a,1);

h_lh = tranproc([h' der],g);
der  = abs(h_lh(:,2));
h_lh = h_lh(:,1);

% Normalization + transformation of t and h 

pos  = 2/(1+alfa);    % inverse of a fraction of positive maxima
cons = 2*pi^4*pos/sqrt(2*pi)/m4/sqrt((1-alfa^2));
                                        %Tm=2*pi*sqrt(m0/m2)/alpha; %mean period between positive maxima

t_lh = t;
h_lh = sqrt(m0)*h_lh;


% Computation of the distribution
[T,H] = meshgrid(t_lh(2:b),h_lh);
f_th  = zeros(a,b);
der   = der(:);
%f_th(:,2:b)=const*der(:,ones(1,b-1)).*(H.^2./(T.^5)).*....
%    exp(-(H./(eps4*T.^2)).^2/8.*((T.^2-alpha^2).^2+....
%    alpha^4*beta^2 ))/(Tm/4*sqrt(m0));
f_th(:,2:b)=cons*der(:,ones(1,b-1)).*(H.^2./(T.^5)).*....
    exp(-0.5*(H./T.^2).^2.*((T.^2-pi^2*m2/m4).^2/(m0*(1-alfa^2))+....
    pi^4/m4));

%pdfplot
f   = createpdf(2);
f.f = f_th;

f.x{1} = t;
f.x{2} = h;
f.labx = {'Tc','Ac'};
f.title='Joint density of (Tc,Ac) - Cavanie et al. (1976)';
%f.eps2=eps2;
f.eps4=eps4;
f.pl = [10:20:90 95 99 99.9];
f.cl = qlevels(f.f,f.pl);
%pdfplot(f)

% Let U,Z be the hight and second deivative (curvature) at a local maximum in a Gaussian proces
% with spectral moments m0,m2,m4. The conditional density ($U>0$) has the following form
%$$
% f(z,u)=c \frac{1}{\sqrt{2\pi}}\frac{1}{\sqrt{m0(1-\alpha^2)}}\exp(-0.5\left(\frac{u-z(m2/m4)}
% {\sqrt{m0(1-\alpha^2)}}\right)^2)\frac{|z|}{m4}\exp(-0.5z^2/m4), \quad z<0, 
%$$
% where $c=2/(1+\alpha)$, $\alpha=m2/\sqrt{m0\cdot m4}$. 
%
% The cavanie approximation is based on the model $X(t)=U \cos(\pi t/T)$, consequently
% we have $U=H$ and by twice differentiation $Z=-U(\pi^2/T)^2\cos(0)$. The variable change has Jacobian
% $2\pi^2 H/T^3$ giving the final formula for the density of $T,H$
%$$
% f(t,h)=c \frac{2\pi^4}{\sqrt{2\pi}}\frac{1}{m4\sqrt{m0(1-\alpha^2)}}\frac{h^2}{t^5}
%       \exp(-0.5\frac{h^2}{t^4}\left(\left(\frac{t^2-\pi^2(m2/m4)}
% {\sqrt{m0(1-\alpha^2)}}\right)^2+\frac{\pi^4}{m4}\right)). 
%$$
%
%
