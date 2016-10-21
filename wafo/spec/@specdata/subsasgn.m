function self = subsasgn(self,index,b,varargin)
% SUBSASGN Define field name assigment for SPECDATA objects
%
% CALL  self = subsasgn(self,index,b,varargin)
%
% Examples
% f = specdata(jonswap);
% f(2) = specdata(torsethaugen);
% [f.data] = deal(3); % replace data with 3
% f.data
% f(2).data= 7;
% f.data
% 
% See also specdata/set

persistent legalFieldName parentFieldName parents

if isempty(legalFieldName) || isempty(parentFieldName) || isempty(parents)
  legalFieldName = fieldnames(self);
  
  parents = parent(self); 
  Np = length(parents);
  parentFieldName = cell(1,Np);
  sself = struct(self);
  for ix=1:Np
    parentFieldName{ix} = fieldnames(sself(1).(parents{ix}));
  end
end


N = length(index);
switch index(1).type
  case '()'
    if N == 1
      if ~isa(b,'specdata')
        error('Not allowed to assign a non-WAFO data object to wdata.')
      end
      self(index(1).subs{:}) = b;
    else
      self(index(1).subs{:}) =  subsasgn(self(index(1).subs{:}),index(2:N),b,varargin{:});
    end
  case '.'
    subs1 = index(1).subs;
    par = cell(1,numel(self));
    [par{:}] = deal(b,varargin{:});
      
    if N >1
      
      switch subs1
        case legalFieldName
          for ix = 1:numel(self)
            self(ix).(subs1) = subsasgn( self(ix).(subs1),index(2:N),par{ix});
          end
        case parentFieldName{1}
           for ix = 1:numel(self)
            self(ix).(parents{1}).(subs1) = subsasgn( self(ix).(parents{1}).(subs1),index(2:N),par{ix});
           end
        case {'S','s'}
          for ix = 1:numel(self)
            self(ix).(parents{1}).data = subsasgn( self(ix).(parents{1}).data,index(2:N),par{ix});
          end
        case {'f','w','theta','k','k2'}
           for ix = 1:numel(self)
             self(ix) = setfwtork(self(ix),index,par{ix});
          end    
      end
    else % N==1
      switch subs1
        case legalFieldName
          if strcmpi(subs1,'wdata') && ~isa(b,'wdata')
            error('Not allowed to assign a non-WAFO data object to wdata.')
          end
          [self.(subs1)] = deal(par{:});
        case parentFieldName{1}
          for ix = 1:numel(self)
            self(ix).(parents{1}).(subs1) = par{ix};
          end
        case {'S','s'}
          for ix = 1:numel(self)
            self(ix).(parents{1}).data = par{ix};
          end
        case {'f','w','theta','k','k2'}
          for ix = 1:numel(self)
             self(ix) = setfwtork(self(ix),index,par{ix});
          end        
        otherwise
        %case {'name','type','data','args','argMin','argMax', 'sampleRate', 'labels',    'title',...
        %    'contourLevels', 'percentLevels',  'workspace',    'info',   'note', 'date'}
        
        try
          self.wdata = subsasgn(self.wdata,index,par{:});
        catch
          error('???? Reference to non-existent field or method %s',subs)
        end
      end
    end
  otherwise
    error('WAFO:SPECDATA:SUBSASGN','???Illegal reference')
end
    
%% subfunction   

 function     self = setfwtork(self,index,b)
   
   N = length(index);
   subs1 = index(1).subs;
   
   x = get([self.wdata],'args');
   is2D = ismember(self.type,{'k2d','dir','encdir'});
   if is2D ~= iscell(x)
     error('args of WDATA object must be a cellarray for 2D data and vector for 1D data!')
   end
   argMin = get(self.wdata,'argMin');
   argMax = get(self.wdata,'argMax');
   
   if is2D
     if strcmpi(self.type(end-2:end),'dir')
       dim = 2*strcmpi(subs1,'theta') + strcmpi(subs1,self.freqtype);
     else
       dim = 2*strcmpi(subs1,'k') + strcmpi(subs1,'k2');
     end
     if dim ==0
       error('???? Reference to non-existent field or method %s',subs1)
     end
     if strcmpi(index(N).type,'()')
       x{dim}(index(N).subs{:}) = b;
     else
       x{dim} = b;
     end
     argMin(dim) = min(b);
     argMax(dim) = max(b);
   else
     if self.freqtype~=subs1
       error('???? Reference to non-existent field or method %s',subs1)
     end
     if strcmpi(index(N).type,'()')
       x(index(N).subs{:}) = b;
     else
       x = b;
     end
     argMin = min(b);
     argMax = max(b);
   end
   self.wdata = set(self.wdata,'args',x,'argMin',argMin,'argMax',argMax);
   