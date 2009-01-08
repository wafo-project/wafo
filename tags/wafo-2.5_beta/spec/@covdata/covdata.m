function  Hs = covdata(data,varargin)
%COVDATA COVDATA class Constructor
%
% CALL hs = covdata(data,args)
%
%

error(nargchk(0,inf,nargin));
switch nargin
 case 0
   wd      = wdata;
   Hs = class(covdatastruct, 'covdata',wd);
 case 1
   switch class(data)
     case 'covdata'
       Hs = data;
     case 'specdata'
       error('Type casting from specdata class not defined yet.')
     case 'struct'      
       Hs    = parseoptions(covdatastruct,data);
       try
         if strcmpi(data.name,'WAFO Covariance Object')
           wd = wdata(data.wdata);
         else
           error('Struct does not appear to be an Spectral density struct!')
         end
       catch
          % Conversion from old covariance struct
       
          ltype = lagtype(data);
          Hs.lagtype = ltype;
          dim = length(ltype);
          
          if dim<1
             error('Struct does not appear to be a Covariance struct!')
          elseif dim==1
            args = data.(ltype);
          elseif dim==2
            args = {data.(ltype(1)),data.(ltype(2))};

          else
            args = {data.(ltype(1)),data.(ltype(2)),data.(ltype(3))};
          end  
          
          labelTxt = {'Lag [m]','Lag [sec]'};
          labels  = labelTxt(1+(ltype=='t'));
          labels{dim+1} = 'ACF';
          
          if data.norm==0
            titleTxt = 'Auto Covariance Function (ACF)';
          else
            titleTxt = 'Auto Correlation Function (ACF)';
          end
          
          
          
          fn = fieldnames(data);
          
          iy = find(strncmpi('R',fn,1) | strcmpi('stdev',fn));
          
          if any(iy)
            Nf = cellfun(@length,fn(iy));
            iy(Nf==1) = [];
            if any(iy)
            for fni = fn(iy).'
              workspace1.(fni{1}) = data.(fni{1});
            end
            else
              workspace1 = [];
            end
          else
            workspace1 = [];
          end          
          wd = wdata(data.R,args);
          set(wd,'note',data.note,'labels',labels,'title',titleTxt,'workspace',workspace1)
       end
       Hs = class(Hs,'covdata',wd);
     case 'wdata'
       Hs = class(covdatastruct, 'covdata',data);
     otherwise
       Hs = covdata;
       Hs = set(Hs,'data',data);
   end
  otherwise
    error('FIX THIS')
    
end

end % function covdata

%% subfunction

function Hs1 = covdatastruct
%COVDATASTRUCT Struct defining legal names and default values for a COVDATA object

Hs1.name = 'WAFO Covariance Object';
Hs1.type = ''; % covtype
Hs1.lagtype = ''; % 'x','y',and/or 't'
Hs1.h    = inf;
Hs1.tr   = [];
Hs1.phi  = 0.;
Hs1.v    = 0.;
Hs1.norm = 0;
end % function covdatastruct
