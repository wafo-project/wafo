function [f, ssinfo]=createwdata(varargin)
% CREATEWDATA Distribution parameter struct constructor
%
% CALL:  f=createwdata(list)
%
%  f = struct with the following fields
%              name: 'WAFO data object'
%              type: 'data'
%              data: []
%            dataCI: []
%              args: {}
%            argMin: []
%            argMax: []
%        sampleRate: NaN
%            labels: {}
%             title: ''
%     contourLevels: []
%     percentLevels: []
%         workspace: []
%              info: []
%    .note         : note string
%    .date         : creation date and time
%
%
%  Examples: 
%  x = linspace(0,10);
%  f = createwdata('args',x,'data',sin(x)) 
%
% See also  wdata

%Tested on: Matlab 7
%History:
% by pab Sep 2007



 f = struct('name','WAFO data object',...
    'type','data',...
    'data',[],...
    'dataCI',[],...
    'args',[],...
    'argMin',[],...
    'argMax',[],...
    'sampleRate',nan,...
    'labels',[],...
    'title','',...
    'caption','',...
    'contourLevels',[],...
    'percentLevels',[],...
    'workspace',[],...
    'child',[],...
    'info',[],...
    'note','',...
    'date',datestr(now));
  f.args = {};
  f.labels = {};
  ssinfo = @setsamplinginfo;
if nargin>0
  f = parseoptions(f,varargin{:});
  if ~isempty(f(1).data)
    f  = setsamplinginfo(f);
  end
end

function  self = setsamplinginfo(self)
%SETSAMPLINGINFO Set argmax, argmin and samplerate
%
% CALL: self = setsamplinginfo(self)

Nf = numel(self);
if Nf>1
  for ix = 1:Nf
    self(ix) = setsamplinginfo(self(ix));
  end
else
try
  if iscell(self.args) % Multidimensional data
    argsz = size(self.args{1});
    isVectorArguments = (sum(argsz>1) <= 1);
    if isVectorArguments
      self.argMax = cellfun(@max,self.args);
      self.argMin = cellfun(@min,self.args);
      dX = cellfun(@diff,self.args,'UniformOutput',false);
      [dX{cellfun(@isempty,dX)}] = deal(nan);
      self.sampleRate = cellfun(@min,dX);
    else
      Nargs       = length(self.args);
      self.argMax        = zeros(1,Nargs);
      self.argMin        = zeros(1,Nargs);
      self.sampleRate  = zeros(1,Nargs);
      if Nargs>3
        % Args assumed made from ndgrid.
        dims = 1:Nargs;
      else
        % Args assumed made from meshgrid.
        dims = [2 1,3:Nargs];
      end
      for ix = 1:Nargs
        self.argMax(ix) = max(self.args{ix}(:));
        self.argMin(ix) = min(self.args{ix}(:));
        dX = diff(self.args{ix},1,dims(ix));
        self.sampleRate(ix) = min(dX(:));
      end
    end    
%     switch length(self.args);
%       case 2,
%         self.type = data_2d;
%       case 3,
%         self.type = data_3d;
%       otherwise 
%         self.type = data_nd;
%     end
  else % One dimensional data
    self.argMax = max(self.args);
    self.argMin = min(self.args);
    self.sampleRate = min(diff(self.args));
    if isempty(self.sampleRate)
      self.sampleRate = nan;
    end
%     self.type = data_1d;
  end
catch
  % Do nothing
end
end

