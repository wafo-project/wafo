function  Hs = fdata(varargin)
%FDATA Distribution parameter class Constructor
%
% CALL hs = fdata(data,args)
%
% FDATA is a container for Distribution parameter objects in WAFO
% 
% Example
%
%
% See also 


% Tested on: matlab 7
% History:
% revised pab updated help header
% By pab sept 2007



thisClass = mfilename;
switch nargin
  case 0
    Hs = class(defaultstruct,thisClass);
  case 1
    switch class(varargin{1})
      case {thisClass}
        Hs = varargin{1};
      case 'struct'
        % Convert from struct to this class
        
        Hs = defaultstruct;
        Hs = parseoptions(Hs(ones(size(varargin{1}))),varargin{1});
        Hs = class(Hs,thisClass);
       
      otherwise
       error('Illegal input!')
    end
  otherwise
    Hs = parseoptions(defaultstruct,varargin{:});
    Hs = class(Hs,thisClass);
end


 function f=defaultstruct
% defaultstruct Struct defining legal names to THIS object
%
% CALL:  f=defaultstruct
%

%Tested on: Matlab 7
%History:
% by pab Sep 2007

f = createfdata;

