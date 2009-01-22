function [ax11, h11, h22  ]=plotmarg2dmom(V,H,phat,condon ,res)
%PLOTMARG2DMOM Plot conditional mean and standard deviation.
%
% CALL: plotmarg2dmom(x1,x2,phat,condon,res);
% 
%       x1,x2 = data
%       phat  =  parameter structure array
%     condon  = 1 mean and standard deviation of X2 given X1
%               2 mean and standard deviation of X1 given X2 (default)
%       res   = resolution (default range(x2)/12)
% 
% 
%  PLOTMARG2DMOM plots the  empirical conditional mean and standard deviation of
%  X1 given X2 or X2 given X1 and compares it with  theoretical quantities.
%
% Example:
%  x1=linspace(0,10)';
%  phat = createfdata('distribution',@pdfmarg2d,'params',[2 2 .5])
%  phat.pdfoptions.distribution={'pdfray','pdfray'};
%  phat.pdfoptions.numpar =ones(1,2);
%  [y1,y2] = rndmarg2d(2000,phat);
%  plotmarg2dmom(y1,y2,phat,2); 
% 
%  See also  mommarg2d


%  tested on: matlab 5.2
% history
% revised pab 03.11.2000
% changed var(tmp) to std(tmp)^2
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vector to structure
% by Per A. Brodtkorb 31.01.99


error(nargchk(3,5,nargin))

if nargin <4 ||isempty(condon), condon =2; end

if nargin<5||isempty(res), 
  if condon==2,  
   res=range(H(:))/12;
  else
   res=range(V(:))/12;
  end 
end


Nmesh=40;
if condon==2,  
   grp=floor(H/res)+1; % dividing the data into groups 
   Ngrp=max(grp);
   cvar=linspace(eps, max(H), Nmesh)';
   [m v ]=findmeancovar(V,grp); 
else
  grp=floor(V/res)+1; % dividing the data into groups 
  Ngrp=max(grp);
  cvar=linspace(eps, max(V), Nmesh)';
  [m v ]=findmeancovar(H,grp); 
end 

h1=linspace(res/2, (Ngrp-0.5)*res, Ngrp)'; 
[M1 V1]= mommarg2d(phat,'condon',condon,'cvar',cvar);
  

if 0,
 [ax1 h11 h22]=plotyy(h1,m,h1,sqrt(v));
 set(ax1(1),'nextplot','add' );
 set(h11, 'LineStyle' , 'x')
 axes(ax1(1));
 plot(cvar,M1,'-'),
 set(ax1(1),'nextplot','replace' );

 set(ax1(2),'nextplot','add' );
 set(h22, 'LineStyle' , 'o')
 axes(ax1(2));
 %axis([0 inf 0 max(v)])

 plot(cvar,sqrt(V1),'-'), 
 set(ax1(2),'nextplot','replace' );
else
  plot(h1,m,'bo',h1,sqrt(v),'rx'); hold on
  plot(cvar,M1,'b-'),
   plot(cvar,sqrt(V1),'r--'),hold off
end
if condon ==2
  xlabel('x2')
else
 xlabel('x1')
end
 title('Conditional mean and standard deviation')
 if nargout>0
   ax11=ax1;
 end


function [m, v ]=findmeancovar(V,grp) 

Ngrp=max(grp);
m=zeros(Ngrp,1);v=m;
for ix=1:Ngrp,
  tmp=V(grp==ix);%find data in group number ix
  
  if length(tmp)>max(4,0),% if less than 4 observations in the group 
    m(ix)=mean(tmp); % mean of data in group ix
    v(ix)=std(tmp).^2;  % variance of data in group ix
  else 
    m(ix)=NaN;
    v(ix)=NaN;
  end
end
return
