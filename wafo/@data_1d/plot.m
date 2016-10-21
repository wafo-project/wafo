function H1=plot(self,f,flag1,flag2,varargin)
%PLOT Plot WDATA objects of DATA_1D type.
%
% CALL: H = plot(self,f,plotflag,sym)
%
%  H = handle to the created object
%
%  plot a struct object with the fields
%      
%      f.data   = matrix
%      f.args   = vector of values 
%  
% optional fields:   
%      f.dataCI
%      f.labels        = cellarray of label strings          (n=1:2) 
%      f.title         = title string
%
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
%      sym      = plot symbol (default '-')  
%
% Example
%  x = linspace(0,10).';
%  f = wdata(sin(x),x);
%  plot(f.type,f,'g')
% 
% See also wdata/plot

% Note: is only able to handle 1D plot i.e. ndim=1

%Tested on: Matlab 5.3, 5.2
%History:
% by pab March 2007
% - based on old pdfplot





hold_state = ishold; % remember old hold state
Nff=length(f);
if Nff>1
 error('WAFO:DATA_1D:PLOT','Illegal size of WDATA object')
end

%cax  = newplot; % axes
%cfig = get(cax,'Parent'); %cfig=gcf;
if iscell(f.args)
  dim = length(f.args);
else
  dim = 1;
end
if dim>length(size(squeeze(f.data)))
  fsiz=size(f.data);
  dim=length(fsiz)-sum(fsiz==1); % Number of non-singleton dimensions
end
if dim~=1
  error('WAFO:DATA_1D:PLOT','Wrong dimension of WDATA object')
end
%dim,
% initialize output to default values
plotflag = 1;
sym='k'; % Black solid line is default

if nargin>2 && ~isempty(flag1)
  if ischar(flag1)
    sym = flag1;
  elseif isnumeric(flag1)
    plotflag = flag1;
  end
end
if nargin>3 && ~isempty(flag2)
  if ischar(flag2)
    sym = flag2;
  elseif isnumeric(flag2)
    plotflag = flag2;
  end
end

if plotflag==0, 
  return,
end
data = transformdata(f.args,f.data,plotflag);
if ~isempty(f.dataCI)
  dataCI = transformdata(f.args,f.dataCI,plotflag);
else
  dataCI = [];
end
H = plot1d(f.args,data,dataCI,sym,plotflag,varargin{:});
    
labelfig(f)
if ~hold_state, 
   hold off, 
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
scale = plotscale(data_1d,plotflag);
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







