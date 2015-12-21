function [value,varargout] = get(self,varargin)
% GET Access data stored in class object
% 
% CALL  values = get(Hs,par1,par2,...)
%
%  Hs         = class object
%  par1,par2..= strings identifying the property to access
%  values     = value of the specified property of the HS object or a
%               struct array containing the specified properties of 
%               the HS object(s).
%
% Examples
%   d = wdata(2,1);
%   get(d)    % displays all properties of d and their values.  
%   val     = get(d,'params')         % data matrix
%   vals    = get(d,'params','distribution')  % struct with data and args
%   allvals = get(d);               % struct with all properties
% %PDFFIT array 
%   d2 = d; d2(2) = wdata(3,2);
%   [val2{1:2}] = get(d2,'data')        % cellarray of data
%   vals2       = get(d2,'data','args') % struct array with data and args
%   allvals2    = get(d2);              % struct array with all properties
%   
%   See also set

%Tested on: matlab 7.3
%History
% By pab 2007
  

n = numel(self);
varargout = cell(1,nargout-1);

switch nargin
  case {1}
    value = struct(self);
  case {2}
    no = max(min(n,nargout),1); 
    [value, varargout{1:no-1}] = deal(self(1:no).(varargin{1}));
  otherwise
    for ix = 1:nargin-1
      [value(1:n).(varargin{ix})] = deal(self.(varargin{ix}));
    end
end
