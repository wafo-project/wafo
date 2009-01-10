function plotweib2dcdf(V,H,varargin) 
%PLOTWEIB2DCDF Plot conditional empirical CDF of X1 given X2=x2
%               and optionally compares it with weib2d distribution.
%  
%  CALL: plotweib2dcdf(x1,x2,a1,b1,a2,b2,c12,options);
%        plotweib2dcdf(x1,x2,phat,options);
% 
%       x1,x2 = data
%   A1 B1,
%   A2 B2 C12 = distribution parameters
%        phat = Distribution parameter struct
%                 as returned from FITWEIB2D.  
%     options = structure with fields:   
%       .condon   :
%       .resolution : resolution (default range(x2)/12)
%       .plotflag : 0  no plotting
%                   1 plot cdf F(x1|x2)                  (default)
%                   2 plot 1-F(x1|x2) on a semilog y-scale 
%                   3 plot F(x1|x2) on a log log scale
%       .rows    : number of rows in each figure    (default 3)
%       .cols    : number of columns in each figure (default 3)
%       .numfigs : total number of figures.         (default 1)
%         sym = {s1,s2} cell array of plot symbols for 
%               plotting empirical and theoretical cdf, respectively.
%               (default {'b.','r--'})
%  
%  PLOTWEIB2DCDF plots the  empirical CDF of X1 given X2 or 
%  X2 given X1 and compares it a 2D Weibull distribution with parameters
%  given by phat.
% 
% NOTE:  SYM can be given anywhere after X1 and X2
% 
% Example:
%  phat = {2 2 2 3 .9};
%  [R1 R2]= rndweib2d(phat{:},1000,1); 
%  plotweib2dcdf(R1,R2,phat{:},'plotflag',2);
%
%  See also  cdfweib2d, plotedf


%  tested on: matlab 5.2
% history
% by pab 25.10.2000


% default values
%~~~~~~~~~~~~~~~
sym = {'b.','r--'}; % default plot symbols for the empirical
                              %  theoretical pdf,respectively

ih = ishold; % save hold state


error(nargchk(2,15,nargin))
Np = 5;
options = struct('condon',2,'resolution',[],'plotflag',1,'rows',3,'cols',3,'numfigs',1); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a1 =params{1};


res = options.resolution;
flag = options.plotflag;

if isempty(res)
  if options.condon==1,
    res=range(V(:))/12;
  else
    res=range(H(:))/12;
  end
end

if flag<1,  return,end

row  = options.rows;
col  = options.cols;
Nfig = options.numfigs;

Nmesh=40;
v1=[];cdfgH=[];
if options.condon==2,
  
  Xc      = V;
  grp     = floor(H/res)+1; % dividing the data into groups 
  Ngrp    = max(grp);
  h1      = linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
  if ~isempty(a1)
    v1      =linspace(eps,max(V)+range(V)/4,Nmesh);
    [X1,X2] = meshgrid(v1,h1);
    cdfgH   = cdfweib2d(X1,X2,params{:},options);
  end
  %max(cdfgH')
  %min(cdfgH')
  xmax    = max(V);
 
else
  Xc      = H;
  grp     = floor(V/res)+1; % dividing the data into groups 
   Ngrp    = max(grp);
  h1      = linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
  if ~isempty(a1)
    v1      = linspace(eps,max(H)+range(H)/4,Nmesh);
    [X1,X2] = meshgrid(v1,h1);
    cdfgH   = cdfweib2d(X2,X1,params{:},options);
  end
  xmax    = max(H);
end


%xmax=min(xmax,4);

fignr    = gcf;
fignrold = fignr;
%figure(fignr)

iy=0;
for ix=Ngrp:-1:min(grp),
  if iy==row*col,iy=0; 
    fignr=fignr+1;
    if Nfig<=fignr-fignrold, break, end
    figure(fignr);
  end
  tmp = Xc(grp==ix);%find data in group number ix
  if length(tmp)>max(3,0),% if less than 6 observations in the group 
    iy=iy+1;
    subplot(row,col,iy)
    if ih, hold on, end % make sure we have hold on for each subplot
    if ~isempty(a1)
      plotedfcnd(tmp,0,[v1; cdfgH(ix,:)]',flag,sym)
    else
      plotedfcnd(tmp,0,[],flag,sym)
    end
    %grid on
    
    axis square
    Ns = 2;
    switch  flag,
      case 1,
	ylabel(['F(x1|x2=' num2str(h1(ix),Ns) ')'])
	xlabel('x1')
	title('')
	axis([0 ceil(xmax) 0 1])
	% title(['Cumulative density function v given h=' num2str(h1(ix),3) ])
      case 2,
	%figtext(0.1,0.1,['1-F(v| h=' num2str(h1(ix),Ns)  ')'],'norm');
	ylabel(['1-F(x1|x2=' num2str(h1(ix),Ns)  ')'])
	xlabel('x1')
	grid off
	title('')
	axis([0 ceil(xmax) 1e-4 1])
	%  title(['The probability of exceeding v given h=' num2str(h1(ix),3)  ])
      case 3 ,
	ylabel(['-log(-log(F(x1|x2=' num2str(h1(ix),Ns)  ')))'])
	xlabel('x1')
	title('')
      case 4,
	ylabel(['(-log(F(x1|x2=' num2str(h1(ix),Ns)  ')))'])
	xlabel('x1')
	title('')
	%  title(['T
    end
    
    %if printflag, print -Pmhlaser ; end   %print -Pmhlaser
  end
end

if ~ih, hold off, end
