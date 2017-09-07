function  Hs = specdata(data,varargin)
%SPECDATA SPECDATA class Constructor
%
% CALL hs = specdata(data,args,type,freqtype,angletype)
%      hs = specdata(wdata,type,freqtype,angletype)    
%
% data      = vector or matrix of spectrum values
% args      = vector or cellarray of vectors of arguments
% wdata     = wdata object
% type      = string defining spectrum type (default 'freq')
% freqtype  = letter defining frequency type (default 'w')
% angletype = string defining angle type of directional spectrum (default 'radians')
% 
% Examples
% Sold = createspec;
% fn   = fieldnames(Sold);
% H = mkjonswap;
% w = linspace(0,4,256)';
% S = specdata(H(w),w);  % Make spectrum object from numerical values
% 
% Sold = get(S,fn{:});    % Convert to old type of spectrum object.
%
% See also wdata covdata

%error(nargchk(0,inf,nargin));
narginchk(0,inf)
switch nargin
 case 0
   Hs = specdatastruct;
   %Hs.data  = wdata;
   Hs = class(Hs, 'specdata',wdata);
 case 1
   switch class(data)
     case 'covdata'
       error('Type casting from covariance object not defined yet')
     case 'specdata'
       Hs = data;
     case 'struct'      
       Hs    = parseoptions(specdatastruct,data);
       try
         if strcmpi(data.name,'WAFO Spectrum Object')
           wd = wdata(data.wdata);
         else
           error('Struct does not appear to be an Spectral density struct!')
         end
       catch
          % Conversion from old spectral density struct
       
          ftype = freqtype(data);
          Hs.freqtype = ftype;
          
          switch lower(Hs.type(end-2:end))
            case 'dir'
              args = {data.(ftype) data.theta};
              if  (abs(max(data.theta)-min(data.theta))>3*pi)
                Hs.angletype = 'degrees';
              else
                Hs.angletype = 'radians';
              end
            case {'req', 'k1d'}
              args = data.(ftype);
            case 'k2d'
              args = {data.(ftype),data.k2};
            otherwise
              error('Struct does not appear to be an old Spectral density struct!')
          end
          
          wd = wdata(data.S,args);
          set(wd,'note',data.note)
       end
       Hs = class(Hs,'specdata',wd);
      try
        Hs = label(Hs);
      catch
      end
      
     case 'wdata'
       Hs = class(specdatastruct, 'specdata',data);
     otherwise
       Hs = class(specdatastruct, 'specdata',wdata(data));
       %Hs = specdata;
       %Hs = set(Hs,'data',data);
   end
  
  case {2,3,4,5}
    Nstart = 0;
    switch class(data)
      case 'wdata'
        Hs = class(specdatastruct, 'specdata',wdata(data));
      case 'double'
        Hs = class(specdatastruct, 'specdata',wdata(data,varargin{1}));
        Nstart = 1;
      otherwise % Numeric values
        error('Unsupported call!')
    end
    ni = nargin-1;
    if ni<Nstart+1 || isempty(varargin{Nstart+1})
      stype = 'freq';
    else
      stype = varargin{Nstart+1};
    end
    if ni<Nstart+2 || isempty(varargin{Nstart+2})
      ftype = 'w';
    else
      ftype = varargin{Nstart+2};
    end
    
    if ni<Nstart+3 || isempty(varargin{Nstart+3})
      if strcmpi(stype(end-2:end),'dir')
        angletype = 'radians';
      else
        angletype = '';
      end
    else
      angletype = varargin{Nstart+3};
    end
    Hs = set(Hs,'type',stype,'freqtype',ftype,'angletype',angletype);
    
  otherwise
    error('FIX THIS')
    
end

end % function specdata

%% subfunction

function Hs1 = specdatastruct
%SPECDATASTRUCT Struct defining legal names and default values for a SPECDATA object

Hs1.name = 'WAFO Spectrum Object';
Hs1.type = ''; % spectype
Hs1.freqtype = '';
Hs1.angletype = '';

Hs1.h    = inf;
Hs1.tr   = [];
Hs1.phi  = 0.;
Hs1.v    = 0.;
Hs1.norm = 0;
end % function specdatastruct
