function [value,varargout] = get(self,varargin)
%SPECDATA/GET Access data stored in a SPECDATA object
% 
% CALL  values = get(Hs,par1,par2,...)
%
%  Hs         = specdata object
%  par1,par2..= strings identifying the property to access
%  values     = value of the specified property of the HS object or a
%               struct array containing the specified properties of the HS object.
%
% Examples
%   d = specdata(jonswap);
%   get(d)    % displays all properties of d and their values.  
%   val     = get(d,'data')          % data matrix
%   vals    = get(d,'data','args')   % struct with data and args
%   allvals = get(d);                % struct with all properties
%   d2 = d; d2(2) = specdata(torsethaugen);
%   [val2{1:2}] = get(d2,'data')        % cellarray of data
%   vals2       = get(d2,'data','args') % struct array with data and args
%   allvals2    = get(d2);              % struct array with all properties
% 
%   See also specdata\set




%Tested on: matlab 7.3
%History
% By pab 2007
  
n = numel(self);
varargout = cell(1,nargout-1);

switch nargin
  case {1}
    value = struct(self);
  case {2}    
   % value = self.(varargin{1});   
   index.type = '.';
   index.subs = varargin{1};
   no = max(min(n,nargout),1); 
   [value,  varargout{1:no-1}] = subsref(self(1:no),index);
   
  otherwise
    
   index.type = '.';
   index.subs = varargin{1};
    for ix = 1:nargin-1
      property = varargin{ix};
      index.subs = property;
      [value(1:n).(property)] = subsref(self,index);
      %value.(property) = self.(property);
    end
end