function [M ,V ,Tm ,cvar, tolm, tolv] = mommarg2d(phat,varargin)
%MOMMARG2D  Mean and variance for the MARG2D distribution.
%
%  CALL:  [M,V] = mommarg2d(phat,options) 
%
%   M,V    = mean and variance of the marg2d distribution
%   phat   = parameter structure array (see fitmarg2d)
%  options = struct with field names
%   condon : 0 returns the mean, covariance and modal value of x1 and x2 (default)
%            1 conditional mean, variance and modal value of x2 given x1 
%            2 conditional mean, variance and modal value of x1 given x2 
%  cvar    : conditional variable i.e., x1 or x2 depending on condon.
%   tol    : relative tolerance (default 1e-3)
%
% Example
%   x1=linspace(0,10)';
%   phat = createfdata('distribution', @pdfmarg2d,'params', [1 2 10]);
%   phat.pdfoptions.distribution = {'pdfray','pdfray'};
%   phat.pdfoptions.numpar = [1 1];
%   [M,V]=mommarg2d(phat,'condon',2,'cvar',x1);
%   plot(x1,M,'r--',x1,sqrt(V),'k-')
%   title(' Conditional mean and standard deviation')
%   legend('E(x1|x2)','std(x1|x2)')
%   xlabel('x2')
% 
% See also  pdfmarg2d, fitmarg2d rndmarg2d

% references:
% Plackett, R. L. (1965) "A class of bivariate distributions."
%                                J. Am. Stat. Assoc. 60. 516-22
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 p. 133-134.


%  tested on: matlab 5.2
% history
  
% revised pab jan2004
%  - added todo comment
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vector to structure
%  Per A. Brodtkorb 28.01.99

% TODO % modal value not implemented yet

options = struct('numpar',[],'distribution',[],'condon',0,'cvar',[],'tol',1e-3);
if nargin==1 && nargout <= 1 && isequal(phat,'defaults')
  M = options;
  return
end
options = parseoptions(options,varargin{:});
tol = options.tol;

pdfoptions = phat.pdfoptions;
try
  VDIST=getdistname(lower(pdfoptions.distribution{1}));
  HDIST=getdistname(lower(pdfoptions.distribution{2}));
catch
   error('Not enough inputs')
end
%cdf_v = str2func(['cdf',VDIST]);
%pdf_v = str2func(['pdf',VDIST]);
mom_v = str2func(['mom',VDIST]);
%cdf_h = str2func(['cdf',HDIST]);
%pdf_h = str2func(['pdf',HDIST]);
mom_h = str2func(['mom',HDIST]);

psi=phat.params(3);
PV=num2cell(phat.params(1:pdfoptions.numpar(1)));
PH=num2cell(phat.params(pdfoptions.numpar(1)+(1:pdfoptions.numpar(2))));


[m1, v1]=mom_v(PV{:});% marginal mean and covariance
[m2, v2]=mom_h(PH{:});% marginal mean and covariance


% modal value not implemented yet
%
 if 0
   % This is a trick to get the html documentation correct.
   k = pdfmarg2d(1,1,2,3);
 end

 switch options.condon
   case 0,    % mean and covariance
     covar=sqrt(v1*v2)*getrho(psi);
     
     M=[m1 m2];
     V=[v1 covar;covar' v2 ];
     %Tm=[tm1,tm2];
     Tm=[];
   case {1,2},
     cvar = options.cvar;
     if options.condon==2,%conditional mean and variance given H
       txt='H';      vr=v2;      mn=m2; 
       %vc=v1;
       mc=m1;
       opt3.condon = 3;
       opt4.condon = 4;
       intfun = @pdfmarg2d;
       
     else %conditional mean and variance given V
       
       txt='V'; vr=v1;     mn=m1;  %vc=v2;
       mc=m2;
       opt3.condon = 13;
       opt4.condon = 14;
       intfun = @localpdfmarg2dswap;
       
     end
     if isempty(cvar),cvar=linspace(0 ,mn+3*sqrt(vr),30)'; end
     if 1
       xinf=mn+mn./mc.*cvar   +15*sqrt(vr); % infinite value for V or H
      %tmp=input(['Enter an infinite value for ', txt, ':  (' , num2str(xinf), ')']);
      %if ~isempty(tmp), xinf=tmp;end
      %disp(['Infinite value for ' ,txt,' is set to ' num2str(xinf(:)')])
     end
   
     M=zeros(size(cvar));V=M;Tm=M;  tolm=M;tolv=V;
     %do the integration with a quadrature formulae
     %E(H|V) or E(V|H)
     [M   tolm]=gaussq(intfun,0,xinf,tol,[],cvar,phat,opt3);
     %E(H^2|V) or E(V^2|H)
     [V  tolv]=gaussq(intfun,0,xinf,tol,[],cvar,phat,opt4);
     V=V-M.^2;%var(H|V) or var(V|H)
     k=find(xinf<M+6*sqrt(V));
     if any(k),
       xinf(k)=M(k)+6*sqrt(V(k))+xinf(k); % update the infinite value
       %disp(['Changed Infinite value for ', txt]);%, ' to ',   num2str(xinf(k))])
        [M(k)   tolm(k)]=gaussq(intfun,0,xinf(k),tol,[],cvar(k),phat,opt3);
        %E(H^2|V) or E(V^2|H)
        [V(k)  tolv(k)]=gaussq(intfun,0,xinf(k),tol,[],cvar(k),phat,opt4);
        V(k)=V(k)-M(k).^2;%var(H|V) or var(V|H)
     end
     
     if nargout>2
       Nint=length(cvar(:));
       
       
       for ix=1:Nint,
         Tm(ix)=fminbnd('-pdfmarg2d(x,P1,P2,P3,P4)',...
           max(0,M(ix)-sqrt(V(ix))),M(ix)+sqrt(V(ix)),[],...
           cvar(ix),phat2,2 ); %modalvalue
       end
       
     end 
 end
end
  function f = localpdfmarg2d(x1,x2,varargin)
  f = pdfmarg2d(x1,x2,varargin{:});
  end
  function f = localpdfmarg2dswap(x1,x2,varargin)
  f = localpdfmarg2d(x2,x1,varargin{:});
  end


% 
% function [me, va]=dist1dstatfun(Ah,dist2 )  
%    switch dist2(1:2)
%       case 'ra',  [me va] = momray(Ah(1));
%       case 'we' ,  [me va]=momweib(Ah(1),Ah(2));
%       case 'gu' ,  [me va]=momgumb(Ah(1),Ah(2),0);
%       case 'tg' ,  [me va]=momgumb(Ah(1),Ah(2),1);
%       case 'ga' ,  [me va]=momgam(Ah(1),Ah(2));
%       case 'lo' ,  [me va]=momlognorm(Ah(1),Ah(2));
%       otherwise, error('unknown distribution')
%     end 
% end

function r=getrho(psi)
r=zeros(size(psi));
k1=find(psi==0);
if any(k1),
 r(k)=-1;
end
k3=find((psi<1.*psi>0)|psi>1);
if any(k3),
 r(k3)=(psi(k3).^2-1-2.*psi(k3).*log(psi(k3)))./((psi(k3)-1).^2) ;
end
end
