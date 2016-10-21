function wafofig5
% WAFOFIG5  Joint distribution (pdf) of crest front velocity and wave height:
%           Theoretical joint density of Vcf and 2*Ac (solid), 
%           kernel density estimate (dash) of Vcf and Hd
%           data from Gullfaks C in the North Sea (dots)
% 
%          ( Increase NNp and Nh to get a smoother theoretical distribution.
%           This may increase the computational time dramatically)

% Adapted to  cssmooth  by GL Feb 2011  
% revised pab Feb2005
% -updated call to kdebin
  
global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global NVcf NHd Nxr Nrate 
global fTcfAc NNp Nh Nnit Nspeed 
global kdeVcfHd Nkernel Nhs NL2

if Nnit<0
  opt = rindoptset('method',abs(Nnit),'speed',Nspeed);
else
  opt = rindoptset('method',0,'nit',(Nnit),'speed',Nspeed);
end
% Only need to calculate Globals which is not empty

% Gullfaks C / North Sea data
% Extract Vcf and Hd
if isempty(NVcf)
  [NVcf NHd] = dat2steep(Nxr,Nrate,0);
end

% Joint distribution of Tcf and Ac
if isempty(fTcfAc)
  Sn=dat2spec(Nxr);
  warning(['This takes several  hours to finish ... (1 hour or more depending' ...
	' on input parameters and your computer)'])
  fTcfAc = spec2thpdf(Sn,0,'TcfAc',[0 11 NNp],Nh,opt);
end
% Transform to the joint distribution of Vcf  and 2*Ac
fVcf2Ac=th2vhpdf(fTcfAc);


% Kernel density estimate of Vcf and Hd
if isempty( kdeVcfHd)
  kopt = kdeoptset('kernel',Nkernel,'hs',Nhs,'L2',NL2);
  kdeVcfHd=kdebin([NVcf, NHd],kopt);
  
  % calculate the levels which encloses fkde.pl percent of the data (v,h)
  r = evalpdf(kdeVcfHd,NVcf, NHd,'linear');
  kdeVcfHd.cl = qlevels2(r,kdeVcfHd.pl);
end  


plot( NVcf, NHd,'.'), hold on
pdfplot( kdeVcfHd,'r--')
pdfplot(fVcf2Ac,'k-')
hold off
wafostamp('Figure 5','(NR)')    

return

function f2 = th2vhpdf(f,v)
% TH2VHPDF Calculates joint density of crest velocity and height 
%         in a stationary Gaussian transform process X(t) where 
%         Y(t) = g(X(t)) (Y zero-mean Gaussian with spectrum given in S). 
%         The transformation, g, can be estimated using lc2tr or dat2tr.
%         
%  CALL:  fvh = th2vhpdf(fth,v);
%
%        fth  = (T Ac) pdf structure 
%        fvh  = (V Ac)=(Ac/T Ac) pdf structure
%          v  = vector of  velocities; note  v >= 0,
%
% Example:
%    fth = spec2thpdf(S,[],'TcfAc',[5 5 51],[],-2,4);
%    fvh = th2vhpdf(fth);  
%
% See also  spec2thpdf

% Tested on : matlab 5.3
% History: 
% revised by Per A. Brodtkorb 19.09.1999

%tic

f2=f;
f2.date=datestr(now);
f2.labx{1}='Velocity (m/s)';
f2.title(find(f2.title=='T'))='V';
k=findstr(f2.title,'min');
if any(k)
 f2.title=[f2.title(1:k(1)) 'ax' f2.title(k(1)+3:end)] ;
end
t=f.x{1}(:);
h=f.x{2}(:);
if nargin<2 || isempty(v)
  v=linspace(0,5,length(t))';
else
  v=v(:);
end
f2.x{1}=v;

for ix=2:length(h)
  v1 =[0; h(ix)./t(2:end) ];
  %hold on,plot(v,f_th.f(ix,2:end).*t(2:end).^2/h(ix) )
  if 0
    % the extrapolation done here does not work well if f.f si sparsely
    % sampled => f2.f values corrupted
    
     v1 =[0; h(ix)./t(2:end) ];
    ind= find(f.f(ix,:)>0);
    f2.f(ix,2:end)=exp(min(cssmooth(v1(ind), log(f.f(ix,ind)'.*t(ind).^2/h(ix)) ,.99,v(2:end) ,1),0));
  else
    % no extrapolation => more robust
     v1 =[ h(ix)./t(2:end) ];
     
    f2.f(ix,2:end)=interp1(v1,[  f.f(ix,2:end)'.*t(2:end).^2/h(ix) ],v(2:end) ,'linear');
  end
end
f2.f(isnan(f2.f))=0;
k=find(f2.f>1000);
if any(k)
  disp('Spurios spikes due to transformation')
  %disp('set them to zero')
  %f2.f(k)=0;
end



switch 2
  case 0, %do nothing
  case 2,% transform Amplitude
    rate=2;
    f2.f=f2.f/rate;
    f2.x{2}=f2.x{2}*rate;
   k=findstr(f2.title,'A');
   if any(k)
     f2.title=[f2.title(1:k(1)-1) num2str(rate) f2.title(k(1):end)] ;
   end 
 case 3
   % this does not work correctly
   g=f.tr;
   utc=f.u;
   
   Ac=f2.x{2}(:);
   Nx=length(Ac);
   Nv=length(f2.x{1}(:));
   der=ones(Nx,1); % dh/dh=1
   hg=tranproc([utc+Ac der],g);
   der1=hg(:,2);
   Acg=hg(:,1); % Gaussian level
   hg=tranproc([-Acg der],fliplr(g));
   der2=hg(:,2);
   der=1+abs(der1.*der2);
   % transform amplitude to wave height H=Ac-G(-g(Ac+utc));
   f2.x{2}=Ac-tranproc(-Acg,fliplr(g)); 
   f2.f=f2.f.*der(:,ones(1,size(f2.f,2)));
   x1=linspace(f2.x{1}(1),f2.x{1}(Nv),Nv)';
   x2=linspace(f2.x{2}(1),f2.x{2}(Nx),Nx)';
   [X1 X2]=meshgrid(x1,x2);
   f2.f = evalpdf(f2,X1,X2,'linear');
   f2.x{2}=x2;
   f2.x{1}=x1;
   
   k=findstr(f2.title,'Ac');
   if any(k)
     f2.title=[f2.title(1:k(1)-1) 'Ac-G(-g(Ac+utc))' f2.title(k(1)+2:end)] ;
   end 
end

if isfield(f,'pl')
  f2.cl=qlevels(f2.f,f.pl);
else
  f2.pl=[10:20:90 95 99 99.9];
  f2.cl=qlevels(f2.f,f2.pl);
end
%toc
