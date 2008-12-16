function [ui, p]=qlevels(pdf,p,x1,x2)
%QLEVELS Calculates quantile levels which encloses P% of PDF
%
%  CALL: [ql PL] = qlevels(pdf,PL,x1,x2);
%
%        ql    = the discrete quantile levels.
%        pdf   = joint point density function matrix or vector
%        PL    = percent level (default [10:20:90 95 99 99.9])
%        x1,x2 = vectors of the spacing of the variables 
%               (Default unit spacing)
%
% QLEVELS numerically integrates PDF by decreasing height and find the 
% quantile levels which  encloses P% of the distribution. If X1 and 
% (or) X2 is unspecified it is assumed that dX1 and dX2 is constant.
% NB! QLEVELS normalizes the integral of PDF to N/(N+0.001) before 
% calculating QL in order to reflect the sampling of PDF is finite.  
% Currently only able to handle 1D and 2D PDF's if dXi is not constant (i=1,2).
%
% Example: 
%   x   = linspace(-8,8,2001);
%   qls = qlevels(pdfnorm(x),[10:20:90 95 99 99.9],x);
% % compared with the exact values
%   ql  = pdfnorm(invnorm((100-[10:20:90 95 99 99.9])/200));  
%   
% See also  qlevels2, tranproc

% Tested on: Matlab 5.3
%History:
% revised pab feb2005
% -updated calls from normpdf to wnormpdf in the help-header example  
% pab 14.10.1999
%  added norm, updated documentation,
% by Per A. Brodtkorb 21.09.1999

 norm=1; % normalize cdf to unity

if any(pdf<0)
  error('This is not a pdf since one or more values of pdf is negative')
end
fsiz=size(pdf);
if min(fsiz)==0,
  ui=[];
  p=[];
  return
end
N=prod(fsiz);
if nargin<3 || isempty(x1) ||((nargin<4) && min(fsiz)>1) || ndims(pdf)>2,
  fdfi=pdf(:);
 
else
  if min(fsiz)==1, % pdf in one dimension
    dx22=1;
  else % pdf in two dimensions
    dx2=diff(x2(:))*0.5;
    dx22=[0 ;dx2]+[dx2;0];
  end
  dx1=diff(x1(:))*0.5;
  dx11=[0 ;dx1]+[dx1;0];
  dx1x2=dx22*(dx11.');
  fdfi=pdf(:).*dx1x2(:);
end


if nargin<2||isempty(p)
  p=[10:20:90 95 99 99.9] ;
elseif any(p<0 | 100<p)
  error('PL must satisfy 0 <= PL <= 100')
end

p2=p/100;

[tmp ind]=sort(pdf(:)); % sort by height of pdf
ind=ind(end:-1:1);
fi=pdf(ind);
fi=fi(:);
%Fi=cumtrapz(fdfi(ind));
Fi=cumsum(fdfi(ind)); % integration in the order of
		      % decreasing height of pdf
		      
if norm  %normalize Fi to make sure int pdf dx1 dx2 approx 1
  Fi=Fi/Fi(end)*N/(N+sqrt(eps));
end
maxFi=max(Fi);
if maxFi>1
  disp('Warning:  this is not a pdf since cdf>1')
  disp('normalizing')
  Fi=Fi/Fi(end)*N/(N+sqrt(eps));
elseif maxFi<.95
  disp('Warning: The given pdf is too sparsely sampled ')
  disp('since cdf<.95.  Thus QL is questionable')
elseif maxFi<0
  error('');
end


ind=find(diff([Fi;1])>0);% make sure Fi is strictly increasing by not considering duplicate values
ui=tranproc(p2(:),[Fi(ind) fi(ind)]); % calculating the inverse of Fi to find the index
                                       % to the desired quantile level
%ui=smooth(Fi(ind),fi(ind),1,p2(:),1) % alternative
%res=ui-ui2


if any(ui>=max(pdf(:)))
  disp('Warning: The lowest percent level is too close to 0%')
end
if any(ui<=min(pdf(:)))
   disp('Warning: The given pdf is too sparsely sampled or')
   disp('         the highest percent level is too close to 100%')   
   ui(ui<0)=0; 
end

return



