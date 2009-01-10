function H1=pdfplot(f,varargin)
%PDFPLOT Plot contents of pdf structures
%
% CALL: H = pdfplot(f,plotflag,x1,x2,x3,sym,method,shading)
%
%  H = handle to the created object
%
%  plot a PDF struct with the pollowing fields
%      
%      f.f      = pdf 
%      f.x      = cellarray of values in n dimensions (n=1:3) 
%  
% optional fields:   
%      f.labx   = cellarray of label strings          (n=1:3) 
%      f.title  = title string
%      f.cl     = contour levels for 2D PDF
%      f.pl     = Percent levels the given contour
%                 lines encloses.  
%     1D:
%      plotflag =  Integer giving plot type, transformation of
%                   data and plot scale defined as: 
%                PlotType  = mod(plotflag,10);
%                TransType = mod(floor(plotflag/10),10);
%                logXscale = mod(floor(plotflag/100),10)>0;
%                logYscale = mod(floor(plotflag/1000),10)>0;
%                logZscale = mod(floor(plotflag/10000),10)>0;  
%          PlotType options are:
%                 1  linear plot of T(x)        (default)
%                 2  stairs plot of T(x)
%                 3  stem plot of T(x)
%                 4  errorbar plot of T(x) (requires nonempty dataCI)%
%                 5  bar plot of T(x)
%                 6  area plot of T(x)
%          where T(x) is the transformed data according to TransType.
%          TransType options are:
%                 0 T(x) = f.data, no transformation
%                 1 T(x) = 1-f.data
%                 2 T(x) = F     = cumtrapz(f.args,f.data)
%                 3 T(x) = 1-F(x) 
%                 4 T(x) = -log(1-F(x)) 
%                 5 T(x) = 10*log10(f.data) 
%     2D:
%      plotflag = 1 contour plot (default)
%                 2 mesh
%                 3 surf
%                 4 waterfall
%                 5 pcolor
%                 6 contour3
%     3D: 
%      plotflag = 1 sliceomatic   (default)
%                 2 sliceomatic (with index x-,y-and z-labels)
%                 3 slice       
%                 4 contour f(X1,X2,X3(x1)), where x1 is an integer
%                 5 contour f(X1(x1),X2,X3), where x1 is an integer
%                 6 contour f(X1,X2(x1),X3), where x1 is an integer
%      x1,x2,x3 = are vectors defining where to slice for 3D data
%                 (default along the axis where f has its maximum)
%      sym      = plot symbol (default '-')  
%      method   = interpolation method for 3D slice 
%                 'linear' (default), 'cubic', or 'nearest'
%      shading  = controls the color shading of SURFACE and PATCH
%                 objects.  'faceted' (default), flat or 'interp'
%
%  Note: - sym,method and shading can be given anywhere after f and in
%          any order.
%
% See also   datastructures, qlevels, cltext

% Note: is only able to handle 1D,2D and 3D plot i.e. ndim=3

%Tested on: Matlab 5.3, 5.2
%History:
% revised pab 
% -removed point-and-click editing because it is obsolete
% revised pab March 2005
% -fixed some bugs
% revised pab 27.11.2002
% -added sliceomatic call -> reordered the plotflag order
% revised pab 07.01.2001
%  - fixed a bug for slice option 3D, reordered plotflag options for 1D.
% revised pab 23.11.2000
% - fixed a bug for calculation of contourlevels when f is not a pdf
% revised pab 08.11.2000
% - added the possibility that f is an array of structs
% - added plotflag 2,3 and 4 for 1D 
% revised pab 15.05.2000
%  - added shading, contour3
%  - now slicing along the axis where f has its maximum
% % revised pab 28.01.2000
%  - added point-and-click editing of all the text objects (title,
%    xlabel, ylabel) of the current figure
%  - improved the printing of contour level text and moved it into a
%    separate function cltext (this may be further improved)
%  - changed see also line
% revised es 24.01.2000 - dim=length(f.x) improving old, added som missing ;   
% revised pab 20.01.2000
%  - added pcolor
% revised pab 18.01.2000
%  - added return statement if plotflag==0
%  - added slice for 3D visualization
%  - added hold_state
%  - added pdfplotchk
% revised pab 5.11.1999
%  - changed PL to pl and CL to cl
% by pab 12.08.99

if ~isstruct(f) % secret option plot a matrix assuming the first column
  % is the independent variable
  f = createpdf('f',f(:,2:end),'x',{f(:,1)});
  
%   [plotflag,sym] = pdfplotchk(varargin,1);
%   switch plotflag
%     case 0, return
%     case 1,  H = plot(f(:,1),f(:,2:end),sym);
%     case 2,  H = semilogy(f(:,1),1-f(:,2:end),sym);
%     case 3,  H = loglog(f(:,1),-log(1-f(:,2:end)),sym);
%     case 4,  H = semilogy(f(:,1),f(:,2:end),sym);
%     case 5,  H = loglog(f(:,1),-log(f(:,2:end)),sym);
%     case 6,  m = size(f,2);
%              H = waterfall(f(:,1),1:m,f(:,2:end));
%     case 11, H = plot(f(:,1),cumtrapz(f(:,1),f(:,2:end)),sym);  
%     case 12, H = semilogy(f(:,1),1-cumtrapz(f(:,1),f(:,2:end)),sym);
%     case 13, H = loglog(f(:,1),-log(1-cumtrapz(f(:,1),f(:,2:end))),sym); 
%     otherwise, error('Unknown option for plotflag')
%   end
%   axis('square')
%   set(gca,'FontSize',12)
%   wafostamp;
%    if (nargout>=1),     H1=H;   end
%   return
end

hold_state = ishold; % remember old hold state
Nff=length(f);
if Nff>1
  cfig=gcf;
  for ix=1:Nff,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    pdfplot(f(ix),varargin{:})
  end
  return
end

cax  = newplot; % axes
cfig = get(cax,'Parent'); %cfig=gcf;

dim = length(f.x);
if dim>length(size(squeeze(f.f)))
  fsiz=size(f.f);
  dim=length(fsiz)-sum(fsiz==1); % Number of non-singleton dimensions
end
%dim,
[plotflag,sym,method,shad,x1,x2,x3] = pdfplotchk(varargin,dim,f);
if plotflag==0, return,end
switch dim			       
  case 1, %1D
    data = transformdata(f.x{1},f.f,plotflag);
    %  if ~isempty(f.dataCI)
    % dataCI = transformdata(f.args,f.dataCI,plotflag);
    %else
    dataCI = [];
    %end
    H = plot1d(f.x{1},data,dataCI,sym,plotflag);
        
  case 2  %2D
    switch plotflag
      case {1,6,7,8,9},
        PL=0;
        if isfield(f,'cl')&&~isempty(f.cl) % check if contour levels is submitted
          CL=f.cl;
          if isfield(f,'pl'),PL=~isempty(f.pl);end % levels defines quantile levels? 0=no 1=yes
        else
          CL=max(f.f(:))-range(f.f(:))*(1-[0.01 0.025 0.05 0.1 0.2 0.4 0.5 0.75]);
          if 0 % automatic levels by using contours
            c=contours(f.x{:},f.f.'); % calculate 8 levels
            if isempty(c)
              c=contours(f.x{:},f.f);%,7); % calculate levels
            end
            %CL = clevels(c);
            limit = size(c,2);
            ix = 1;
            while(ix < limit)
              CL(ix) = c(1,ix);
              npoints = c(2,ix);
              nexti = ix+npoints+1;
              c(:,ix)=NaN;
              ix = nexti;
            end
            CL=unique(CL);
          end
        end
        if PL,
          clvec=sort(f.pl);
        else
          clvec=sort(CL);
        end
        if any(plotflag==[1 8 9])
          [cs hcs] = contour(f.x{:},f.f,CL,sym);
        else
          [cs hcs] = contour3(f.x{:},f.f,CL,sym);
        end
        if any(plotflag==[1,6])
          ncl=length(clvec);
          if ncl>12, ncl=12; disp('   Only the first 12 levels will be listed in table.'),end
           axcl = cltext(clvec(1:ncl),PL);  % print contour level text
        elseif any(plotflag==[7 9])
          clabel(cs);
        else
          clabel(cs,hcs);
        end
	
      case 2,	mesh(f.x{:},f.f); % meshz
      case 3,	surf(f.x{:},f.f);  %shading interp % flat, faceted       % surfc
      case 4,	waterfall(f.x{:},f.f);
      case 5, pcolor(f.x{:},f.f); %shading interp % flat, faceted
      case 10,
        [cs,hcs]=contourf(f.x{:},f.f); clabel(cs,hcs); fcolorbar(cs);
      otherwise, error('unknown option for plotflag')
    end
    if any(plotflag==(2:5))
       shading(shad);
    end
  case 3, %3D
    switch plotflag
     case 1,
      sliceomatic(f.x{:},f.f)
     case 2,
      sliceomatic(f.f)
     case 3,
      [X,Y,Z]=meshgrid(f.x{:});
      %method='linear';%, 'cubic','spline', or 'nearest'
      slice(X,Y,Z,f.f,x1,x2,x3,method);
      shading(shad);
     case 4,
      x1 = round(x1);
      contour(f.x{1:2},f.f(:,:,x1));
      if isempty(f.labx{3}),f.labx{3}='x3';end
      f.title=[f.title f.labx{3} ' = ' num2str(f.x{3}(x1))];
      f.labx{3}=[];
     case 5,
      x1 = round(x1);
      contour(f.x{[2 3]},squeeze(f.f(:,x1,:)).');
      if isempty(f.labx{1}),f.labx{1}='x1';end
      f.title=[f.title f.labx{1} ' = ' num2str(f.x{1}(x1))];
      f.labx{1}=[];
      f.labx=f.labx([2 3 1]);
     case 6,
      x1 = round(x1);
      contour(f.x{[1 3]},squeeze(f.f(x1,:,:)).');
      if isempty(f.labx{1}),f.labx{2}='x2';end
      f.title=[f.title f.labx{2} ' = ' num2str(f.x{2}(x1))];
      f.labx{2}=[];
      f.labx=f.labx([1 3 2]);
    end
end
if isfield(f,'labx') && max(size(f.labx))>=1
  Nf=max(size(f.labx));
  if Nf>=3, zlabel(f.labx{3}), end
  if Nf>=2, ylabel(f.labx{2}), end
  xlabel(f.labx{1})
end
if isfield(f,'title')
  title(f.title)
end
axis('square')
set(gca,'FontSize',12)
wafostamp;



if ~hold_state, 
   hold off, 
   %set(cfig,'NextPlot', 'replace'); 
end % reset to old hold state

if (nargout>=1)
  H1=H;
end
return

function H = plot1d(args,data,dataCI,sym,plotflag,varargin)

plottype = mod(plotflag,10);
switch plottype
  case 0 %  No plotting
    H = [];
    return
  case 1, H = plot(args,data,sym,varargin{:});
  case 2, H = stairs(args,data,sym,varargin{:});
  case 3, H = stem(args,data,sym,varargin{:});
  case 4, 
    if isempty(dataCI)
      error('WAFO:DATA_1D:PLOT','Must have non-empty dataCI in order to plot errorbars!')
    end
    H = errorbar(args,data,dataCI(:,1)-data(:),dataCI(:,2)-data(:),sym,varargin{:});
  case 5, H = bar(args,data,sym,varargin{:});
  case 6, 
    level = str2double(sym);
    if isfinite(level)
      H = area(args,data,level,varargin{:});
    else
      H = area(args,data,varargin{:});
    end
  otherwise
    error('Unknown plotflag')
end
scale = plotscale(plotflag);
logXscale = any(scale=='x');
logYscale = any(scale=='y');
logZscale = any(scale=='z');
%logXscale = mod(floor(plotflag/100),10)>0;
%logYscale = mod(floor(plotflag/1000),10)>0;
%logZscale = mod(floor(plotflag/10000),10)>0;
if logXscale, set(gca,'Xscale','log');end
if logYscale, set(gca,'Yscale','log');end
if logZscale, set(gca,'Zscale','log');end



% Should probably move this to specdata/plot instead
%fmin = min(data(:));
 transFlag = mod(floor(plotflag/10),10);
 logScale = logXscale || logYscale || logZscale;
if  logScale || (transFlag ==5 && ~logScale)
  ax = axis;
  fmax1 = max(data(:));
  if transFlag ==5 && ~logScale
    ax(4) = 11*log10(fmax1);
    ax(3) = ax(4)-40;
  else
    ax(4) = 1.15*fmax1;
    ax(3) = ax(4)*1e-4;
  end
  axis(ax)
end
if ~isempty(dataCI) && plottype < 3
  hold('on')
  %H2 = plot1d(args,dataCI,[],'r--',plotflag);
  plot1d(args,dataCI,[],'r--',plotflag);
end

function scale=plotscale(plotflag)
% DATA_1D/PLOTSCALE Return plotscale from plotflag
%
% CALL scale = plotscale(plotflag)
%
% plotflag = integer defining plotscale.
%   Let scaleId = floor(plotflag/100). 
%   If scaleId < 8 then:
%      0 'linear' : Linear scale on all axes.
%      1 'xlog'   : Log scale on x-axis.
%      2 'ylog'   : Log scale on y-axis.
%      3 'xylog'  : Log scale on xy-axis.
%      4 'zlog'   : Log scale on z-axis.
%      5 'xzlog'  : Log scale on xz-axis.
%      6 'yzlog'  : Log scale on yz-axis.
%      7 'xyzlog' : Log scale on xyz-axis.
%  otherwise
%   if (mod(scaleId,10)>0)            : Log scale on x-axis.
%   if (mod(floor(scaleId/10),10)>0)  : Log scale on y-axis.
%   if (mod(floor(scaleId/100),10)>0) : Log scale on z-axis.
%
% scale    = string defining plotscale valid options are:
%       'linear', 'xlog', 'ylog', 'xylog', 'zlog', 'xzlog',
%       'yzlog',  'xyzlog' 
%
% Example
% plotscale(data_1d,100)  % xlog
% plotscale(data_1d,200)  % ylog
% plotscale(data_1d,1000) % ylog
%
% See also data_1d/plotscaleflag
scaleId = floor(plotflag/100);
if scaleId<8
  scaleId = scaleId+1;
else
  logXscaleId = (mod(scaleId,10)>0);
  logYscaleId = (mod(floor(scaleId/10),10)>0)*2;
  logZscaleId = (mod(floor(scaleId/100),10)>0)*4;
  scaleId = logYscaleId +logXscaleId+logZscaleId +1;
end
scales = {'linear','xlog','ylog','xylog','zlog','xzlog','yzlog','xyzlog'};

scale = scales{scaleId};

function data = transformdata(x,f,plotflag)
transFlag = mod(floor(plotflag/10),10);
switch transFlag
  case 0, data = f;
  case 1, data = 1-f;
  case 2, data = cumtrapz(x,f);
  case 3, data = 1-cumtrapz(x,f);
  case {4,5},
    msgid = 'MATLAB:log:logOfZero';
    warnstate = warning('query', msgid);
    warning('off', msgid)
    if transFlag==4
      %data = -log(1-cumtrapz(x,f));
      data = -log1p(-cumtrapz(x,f));
    else
      if any(f(:)<0)
        error('Invalid plotflag: Data or dataCI is negative, but must be positive')
      end
      data = 10*log10(f);
    end
    warning(warnstate);
  otherwise
    error('Unknown plotflag')
end





function [plotflag,sym,method,shad,x1,x2,x3] = pdfplotchk(P,dim,f)
%pdfplotCHK Helper function for pdfplot.
%
% CALL  [plotflag,sym,method,shad,x1,x2,x3]=pdfplotchk(P,dim) 
%
%   P = the cell array P of input arguments (between 0 and 6 elements)


% initialize output to default values
plotflag = 1;
sym='k-'; % Black dots is default
method='linear'; % linear is default
shad ='faceted';
if dim==3
  % Old call 
  %x1=mean(f.x{1});x2=mean(f.x{2});x3=mean(f.x{3});  
  % New call slicing where the maximum value is located
  [fmax, ind] = max(f.f(:));
  [I2,I1,I3] = ind2sub(size(f.f),ind);
  x1=f.x{1}(I1);x2=f.x{2}(I2);x3=f.x{3}(I3);
else
  x1=[];x2=[];x3=[];
end

Np=length(P);
try
  strix = cellfun('isclass',P,'char');
catch
  strix=zeros(1,Np);
  for ix=1:Np, % finding symbol strings 
    strix(ix)=ischar(P{ix});
  end
end
k=find(strix);
if any(k) % remove strings
  Nk=length(k);
  if Nk>3
    warning('WAFO:PDFPLOT','More than 3 strings are not allowed in ')
  end
  for ix=1:length(k)
    switch lower(P{k(ix)})
    case {'flat','faceted','interp'}  , shad   = P{k(ix)};
    case {'linear','cubic', 'nearest'}, method = P{k(ix)};
    otherwise % plotsymbol is given
      sym = P{k(ix)};
    end
  end
  Np=Np-length(k);
  P={P{find(~strix)}}; % remove strings from input
end


if (Np>0) && ~isempty(P{1})
  plotflag=P{1};
end
if dim==3
  switch plotflag
    case {1,2,3}, % do nothing
    case 4, x1 = I3;
    case 5, x1 = I1;
    case 6, x1 = I2;
  end
end


if (Np>1) && ~isempty(P{2})
  x1=P{2};
end
if (Np>2) && ~isempty(P{3})
  x2=P{3};
end
 
if (Np>3) && ~isempty(P{4})
  x3=P{4};
end
return