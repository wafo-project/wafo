function H1=plot(self,f,varargin)
%PLOT Plot WDATA objects of DATA_3D type.
%
% CALL: H = plot(self,f,plotflag,x1,x2,x3,sym,method,shading)
%
%  H = handle to the created object
%
%  plot a WDATA object with the fields
%      
%      f.data   = matrix
%      f.args   = vector of values 
%  
% optional fields:   
%      f.dataCI
%      f.labels        = cellarray of label strings          (n=1:2) 
%      f.title         = title string
%
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
% Example
%  x = linspace(0,10).';
%  [X1,X2] = meshgrid(x,x);
%  f = wdata(peaks(X1,X2),{x,x});
%  plot(f.type,f)
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
 error('WAFO:DATA_2D:PLOT','Illegal size of WDATA object')
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
if dim~=3
  error('WAFO:DATA_3D:PLOT','Wrong dimension of WDATA object')
end

[plotflag,sym,method,shad,x1,x2,x3] = plotchk(varargin,dim,f);


if plotflag==0, 
  return,
end

switch plotflag
 case 1,
   if sum(size(f.args{1})>1)>1
     sliceomatic(f.data,f.args{1}(1,:,1),f.args{2}(:,1,1),f.args{3}(1,1,:))
   else
     sliceomatic(f.data,f.args{:})
   end
 case 2,
  sliceomatic(f.data)
 case 3,
   H=slice(f,x1,x2,x3,method);
  shading(shad);
 case 4,
  x1 = round(x1);
  [cs,H] = contour(f.args{1:2},f.data(:,:,x1));
  if isempty(f.labx{3}),f.labx{3}='x3';end
  f.title=[f.title f.labx{3} ' = ' num2str(f.args{3}(x1))];
  f.labels{3}=[];
 case 5,
  x1 = round(x1);
  [cs,H] = contour(f.args{[2 3]},squeeze(f.data(:,x1,:)).');
  if isempty(f.labx{1}),f.labx{1}='x1';end
  f.title=[f.title f.labx{1} ' = ' num2str(f.args{1}(x1))];
  f.labels{1}=[];
  f.labels = f.labels([2 3 1]);
 case 6,
  x1 = round(x1);
  [cs,H] = contour(f.args{[1 3]},squeeze(f.data(x1,:,:)).');
  if isempty(f.labx{1}),f.labx{2}='x2';end
  f.title=[f.title f.labx{2} ' = ' num2str(f.args{2}(x1))];
  f.labels{2}=[];
  f.labels=f.labels([1 3 2]);
end
labelfig(f)
if ~hold_state
  hold('off')
end
if (nargout>=1)
  H1=H;
end
end






function [plotflag,sym,method,shad,x1,x2,x3] = plotchk(P,dim,f)
%plotCHK Helper function for plot.
%
% CALL  [plotflag,sym,method,shad,x1,x2,x3]=plotchk(P,dim) 
%
%   P = the cell array P of input arguments (between 0 and 6 elements)


% initialize output to default values
plotflag = 1;
sym='k-'; % Black dots is default
method='linear'; % linear is default
shad ='faceted';
if dim==3
  % Old call 
  %x1=mean(f.args{1});x2=mean(f.args{2});x3=mean(f.args{3});  
  % New call slicing where the maximum value is located
  [fmax, ind] = max(f.data(:));
  [I2,I1,I3] = ind2sub(size(f.data),ind);
  x1 = f.args{1}(I1);
  x2 = f.args{2}(I2);
  x3 = f.args{3}(I3);
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
    warning('WAFO:WDATA:PLOT','More than 3 strings are not allowed in ')
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
end
