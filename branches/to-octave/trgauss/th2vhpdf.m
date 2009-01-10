function f2 = th2vhpdf(f,v)
%TH2VHPDF Transform joint T-H density to V-H density
%         
%  CALL:  fvh = th2vhpdf(fth,v);
%
%        fth  = (T Ac) pdf structure 
%        fvh  = (V Ac)=(Ac/T Ac) pdf structure
%          v  = vector of  velocities; note  v >= 0,
%
% Example:
%    opt = rindoptset('speed',5,'nit',1,'method',0);
%    S   = jonswap;  
%    fth = spec2thpdf(S,4,'TcfAc',[0 5 51],[],opt);
%    fvh = th2vhpdf(fth);  
%
% See also  spec2thpdf

% Tested on : matlab 5.3
% History: 
% revised pab 3Dec2003  
% revised by Per A. Brodtkorb 19.09.1999

%tic

f2=f;
f2.date=datestr(now);
f2.labx{1}='Velocity (m/s)';
f2.title((f2.title=='T'))='V';
k=findstr(f2.title,'min');
if any(k)
 f2.title=[f2.title(1:k(1)) 'ax' f2.title(k(1)+3:end)] ;
end
t=f.x{1}(:);
h=f.x{2}(:);
if nargin<2||isempty(v)
  v=linspace(0,4,length(t))';
else
  v=v(:);
end
f2.x{1}=v;
f2.f = zeros(length(f2.x{2}),length(v));
k0 = find(v);
k  = find(t);
k3 = find(h(:)).';
method = 'linear'; 
for ix=k3
  if 0
    % the extrapolation done here does not work well if f.f si sparsely
    % sampled => f2.f values corrupted
    
     v1 = h(ix)./t(k) ;
    ind = find(f.f(ix,k)>0);
    f2.f(ix,k0)=exp(min(smooth(v1(ind), log(f.f(ix,k(ind)).'.*t(k(ind)).^2/h(ix)) ,.99,v(k0) ,1),0));
  else
    % no extrapolation => more robust
     v1 = h(ix)./t(k) ;
    f2.f(ix,k0)=interp1(v1,f.f(ix,k).'.*t(k)./v1 ,v(k0) ,method);
  end
end
f2.f(isnan(f2.f))=0;
k=find(f2.f>1000);
if any(k)
  disp('Spurios spikes due to transformation')
  %disp('set them to zero')
  %f2.f(k)=0;
end



switch 0
  case 0, %do nothing
  case 2,% transform Amplitude
    rate=1.9;
    f2.f=f2.f/rate;
    f2.x{2}=f2.x{2}*rate;
   k=findstr(f2.title,'A');
   if any(k)
     f2.title=[f2.title(1:k(1)-1) num2str(rate) f2.title(k(1):end)] ;
   end 
 case 3
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