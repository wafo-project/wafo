function f=trraylpdf(x,def,gn,utc)
%TRRAYLPDF Transformed Rayleigh approximation for amplitudes 
%         
%  CALL:  f = trraylpdf(x,def,g,u);
%
%    f    = density structure of wave amplitude f(x);
%
%    x    = a row vector with x-values.
%    def  = 'Ac',    gives wave crest amplitude Ac (default).
%           'At',    gives wave trough amplitude At.
%           'AcAt',  gives wave range Ac+At
%    g    = [y g(y)] a two column matrix with the transformation  g(y).
%    u    = reference level (default the most frequently crossed level).
%
% Example:
%  np =10000; dt = .2; method = 1
%  x1 = spec2nlsdat(jonswap,np,dt);
%  [S, H,Ac,At] = dat2steep(x1,4,method); 
%  gnl = dat2tr(x1,'nonlinear');   % nonlinear transformation
%  gl  = dat2tr(x1,'linear');      % identity transformations  
%  x   = linspace(0,8);
%  fnl = trraylpdf(x,'ac',gnl,0);  % with transformation
%  fl  = trraylpdf(x,'ac',gl,0);   % without transformation
%  plotedf(Ac), hold on
%  pdfplot(fnl,21,'g'), 
%  pdfplot(fl,21,'r'),  hold off  
%
% See also  pdfray

% Tested on : matlab 5.3
% svi    10.11.1999
% ir     24.05.2000
% ir     25.06.2000
% Revised pab 
% Added example  


% References:
% Rychlik, I. and Leadbetter, M.R. (1997)
% 'Analysis of ocean waves by crossing- and oscillation-intensities'
% Proceedings of Seventh (1997) ISOPE Conference, Vol. III, pp. 206-213.
%

x = x(:);  
f=createpdf;

if nargin<2||isempty(def)
  def='ac';
end

switch lower(def(2:end))
 case  'c',               defnr = 1;
 case  't',               defnr =-1;
 case  'cat',             defnr = 2;
 otherwise, error('Unknown def')
end

if  nargin<3||isempty(gn)
  gn=[(-5:0.02:5)' (-5:0.02:5)'];
end

if nargin<4||isempty(utc)
  utc_d = gaus2dat([0 0],gn); % most frequent crossed level 
  utc   = utc_d(1,2);
end

if  nargin<1||isempty(x)
   if defnr<2   
%      xx=gaus2dat([ (0:0.1:5)' (0:0.1:5)'],gn);
%      x=xx(:,2)-utc;
       ac=gaus2dat([0 5],gn);
       x=(linspace(utc,ac(1,2),100)-utc)';
   else
%      xx1=gaus2dat([ (0:0.05:5)' (0:0.05:5)'],gn);
%      xx2=gaus2dat([ (0:0.05:5)' -(0:0.05:5)'],gn);
%      x=xx1(:,2)-xx2(:,2);
       ac=gaus2dat([0 5],gn);
       at=gaus2dat([0 -5],gn);
       x=linspace(0,ac(1,2)-at(1,2),100)';
   end
end

h1=x;
f.x={x};

der1=ones(length(h1),1);der2=ones(length(h1),1);

if defnr==1
   hg1=tranproc([utc+h1 der1],gn);der1=abs(hg1(:,2));
%   f.f  = der1.*pdfray(hg1(:,1),1);
   f.f  = der1.*hg1(:,1).*exp(-0.5*hg1(:,1).^2);
end

if defnr==-1
  hg2=tranproc([utc-h1 der2],gn);der2=abs(hg2(:,2));
%  f.f  = der2.*pdfray(-hg2(:,1),1);
   f.f  = -der2.*hg2(:,1).*exp(-0.5*hg2(:,1).^2);
end 

if defnr==2
  r=0:0.01:7;
  h11=gaus2dat([ones(length(r),1),r'],gn);%G(R)
  h22=gaus2dat([ones(length(r),1),-r'],gn);%G(-R)
                                             %new g transformation G(R)-G(-R)
  h=h11(:,2)-h22(:,2);
  Gs=[r' h];gs=fliplr(Gs);
                                             %derivative of g
  der=ones(length(r),1);
  hgs=tranproc([h der],gs);
  der=abs(hgs(:,2));  r2=abs(hgs(:,1));
%  ff=der.*pdfray(r2,1);
  ff  = der.*r2.*exp(-0.5*r2.^2);
  f.f  = interp1(h,ff,x,'linear');
end

switch lower(def)
 case  'ac'
  Htxt = 'Transformed Rayleigh approx. of Ac density';
  xtxt = 'Ac [m]';
 case  'at'
  Htxt = 'Transformed Rayleigh approx. of At density';
  xtxt = 'At [m]';
 case  'acat'
  Htxt = 'Transformed Rayleigh approx. of H=Ac+At density';
  xtxt = 'H=Ac+At [m]';
end 
  
f.title=Htxt;
f.labx{1}=xtxt;


