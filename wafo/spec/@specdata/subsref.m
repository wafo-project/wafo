function [b,varargout]= subsref(self,index)
% SUBSREF Define field name indexing or methods to apply on SPECDATA objects 
%
%  CALL [b,varargout]= subsref(self,index)
%
% Examples
%  S = specdata(jonswap);
%  ft = S.freqtype;
%  plot(S.(ft),S.S)            % Plot spectrum
%  S.plot                      % Alternative call
%  argMax = S.argMax;          % = S.wdata.argMax
%  argMax = S.wdata.argMax;
%
%  S2 = S; S2(2) = specdata(torsethaugen);
%  S2.plot       % Plot both spectra in separate figures
%
%  am = {S2.argMax}            % cellarray of argMax values
%  am1 = {S2(1:2).argMax}      % NB!: Matlabbug returns only the first value
%  [am2{1:2}] = S2(1:2).argMax % Use this call instead!
%  
% See also specdata/get


% TODO % This is slow. Optimize this for speed.
% TODO % am1 = {S2(1:2).argMax} is a Matlabbug returns only the first value

%History
%By pab 2007

persistent legalMethodName legalFieldName parentFieldName outPutMethods
if isempty(legalMethodName) || isempty(legalFieldName) || ...
    isempty(parentFieldName) || isempty(outPutMethods)
  if isoctave
    legalMethodName = methods(self);
  else  
    legalMethodName = methods(self, '-full');
  end
  legalFieldName = fieldnames(self);
  
  parents = parent(self);
  %sself = struct(self);
  parentFieldName = cell(0,1);
  for ix=1:length(parents)
    parentFieldName = cat(1,parentFieldName,fieldnames(self(1).(parents{ix})));
  end
  
  outPutMethods = {'char','fieldnames'};
end


N = length(index);

varargout = cell(1,nargout);

subs1 = index(1).subs;
switch index(1).type
  case '()'
    % Working only on a part of self
    pself  = self(subs1{:}); % copy part of self
    if N == 1
      b =  pself;
    elseif strcmpi(index(2).type,'.')
      switch index(2).subs
        case legalMethodName
          N3 = min(N,3);
          if nargout<1 && ~any(strcmpi(index(2).subs,outPutMethods))
            feval(index(2).subs,pself,index(3:N3).subs{:});
            if strncmpi(index(2).subs,'set',3)
              self(subs1{:}) = pself;
              assignin('caller',inputname(1), self);
            end
          else
            nout = min(nargout(index(2).subs),nargout);
            [b,varargout{1:nout-1}] = feval(index(2).subs,pself,index(3:N3).subs{:});
          end
        otherwise
          nout = max(min(numel(pself),nargout),1);
          [b, varargout{1:nout-1}] = subsref(pself,index(2:end));
      end % switch
    else
      error('WAFO:WDATA:SUBSREF','???? Reference to field or method expected.')
    end % if
    
    
  case '{}'
    error('WAFO:WDATA','??? Cell contents reference from a non-cell array object.')
    
    
  case '.'
    Nself = numel(self);
    nout = max(min(Nself,nargout),1);
    switch subs1
      case legalFieldName
       
        if N == 1
          [b,varargout{1:nout-1}] =  deal(self(1:nout).(subs1));
        else
          b = subsref(self(1).(subs1),index(2:end));
          for ix=2:nout
            varargout{ix-1} =  subsref(self(ix).(subs1),index(2:end));
          end
          if nout<2 &&  nargout>1
            [varargout{:}] = deal(b);
          end
        end
        
      case {'S','s'}
        if N==1
          pwdata = struct([self.wdata]);
          [b,varargout{1:nout-1}] =  deal(pwdata.data);
        else
         b = subsref(self(1).wdata.data,index(2:end));          
         for ix=2:nout
           varargout{ix-1} =  subsref(self(ix).wdata.data,index(2:end));
         end
        end
           
      case {'f','w','theta','k','k2'}
        
       b = getfwtork(self(1),index);
       for ix=2:nout
         varargout{ix-1} = getfwtork(self(ix),index);
       end
        
      case parentFieldName
        if N ==1
           pwdata = struct([self.wdata]);
          [b,varargout{1:nout-1}] =  deal(pwdata(1:nout).(subs1));
        else
          b = subsref(self(1).wdata.(subs1).wdata,index(2:end));
          for ix=2:nout
            varargout{ix-1} =  subsref(self(ix).wdata.(subs1),index(2:end));
          end
        end
        
      case legalMethodName
        N2 = min(N,2);
        if nargout<1 && ~any(strcmpi(subs1,outPutMethods))
          feval(subs1,self,index(2:N2).subs{:});
          if strncmpi(subs1,'set',3)
            assignin('caller',inputname(1), self);
          end
        else
          nout = min(nargout(subs1),nargout);
          [b,varargout{1:nout-1}] = feval(subs1,self,index(2:N2).subs{:});
         
        end
      otherwise
        error('WAFO:WDATA','???? Reference to non-existent field or method %s',index(1).subs)
    end % switch
end % switch


 function b = getfwtork(self,index)
        
   N = length(index);
   subs1 = index(1).subs;
   x = get(self.wdata,'args');
   is2D = any(strcmpi(self.type,{'k2d','dir','encdir'}));
   if is2D ~= iscell(x)
     error('args of WDATA object must be a cellarray for 2D data and vector for 1D data!')
   end
   if is2D
     if strcmpi(self.type(end-2:end),'dir')
       dim = 2*strcmpi(subs1,'theta') + strcmpi(subs1,self.freqtype);
     else
       dim = 2*strcmpi(subs1,'k') + strcmpi(subs1,'k2');
     end
     if dim ==0
       error('???? Reference to non-existent field or method %s',subs1)
     end
     b = x{dim};
   else
     if ~strcmpi(self.freqtype,subs1)
       error('???? Reference to non-existent field or method %s',subs1)
     end
     b = x;
   end
   if N>1
     b = subsref(b,index(2:N));
   end
   
return 





%% old call kept just in case!
% b = self;
% for ix = 1:N
%   subs = index(ix).subs;
%   switch index(ix).type
%     case '()'
%       b = b(subs{:});
%     case '.'   
%       if any(strcmp(subs,legalFieldName))
%         if ix == N
%           [b,varargout{1:nargout-1}] = deal(b.(subs));
%           return
%         else          
%           b = b.(subs);
%         end
%       elseif isa(b,'specdata')
%         switch subs
%           case parentFieldName
%             [b,varargout{1:nargout-1}] = get([b.wdata],subs);
%           case {'getspectrum', 's','S'}
%             [b,varargout{1:nargout-1}] = get([b.wdata],'data');
%           case {'f','w','theta','k','k2'}
%             x = get(b.wdata,'args');
%             is2D = any(strcmpi(self.type,{'k2d','dir','encdir'}));
%             if is2D ~= iscell(x)
%               error('args of WDATA object must be a cellarray for 2D data and vector for 1D data!')
%             end
%             if is2D
%               if strcmpi(b.type(end-2:end),'dir')
%                 dim = 2*strcmpi(subs,'theta') + strcmpi(subs,b.freqtype);
%               else
%                 dim = 2*strcmpi(subs,'k') + strcmpi(subs,'k2');
%               end
%               if dim ==0
%                  error('???? Reference to non-existent field or method %s',subs)
%               end
%               b = x{dim};
%             else
%               if ~strcmpi(b.freqtype,subs)
%                 error('???? Reference to non-existent field or method %s',subs)
%               end
%               b = x;
%             end
%             
%           otherwise
%             
%             if (ix == 1) || (ix==2)
%               if (N==2 || N==3) && strcmpi(index(N).type,'()')
%                 
%                 iz = find(strcmpi(subs,legalMethodName));
%                 if any(iz)
%                   if nargout<1 && ~any(strcmpi(legalMethodName{iz},outPutMethods))
%                     feval(legalMethodName{iz},b,index(ix+1).subs{:});
%                     if isa(b,'specdata')
%                       % Write the newly updated data object back to the calling program
%                       assignin('caller',inputname(1), b);
%                     elseif isa(b,'wdata')
%                        self.wdata = b;
%                       assignin('caller',inputname(1), self);
%                     end
%                     clear b
%                   else
%                     [b,varargout{1:nargout-1}]  = feval(legalMethodName{iz},b,index(ix+1).subs{:});
%                   end
%                 else
%                   error('???? Reference to non-existent field or method %s',subs)
%                 end
%                 return
%               else
%                 error('???? Reference to non-existent field %s',subs)
%               end 
%             
%             end
%         end % switch
%         
%       else
%         try
%           if nargout>1
%               [b varargout{1:nargout-1}] = subsref(b,index(ix:N));
%           elseif nargout<1 && strcmpi(index(N).type,'()') && ...
%                  any(strcmpi(index(N-1).subs,legalMethodName))
%              subsref(b,index(ix:N))
%              clear b  
%           else            
%             b = subsref(b,index(ix:N));
%           end
%           
%           return
%         catch
%           error('???? Reference to non-existent field or method %s',subs)
%         end
%       end % if
%   end
% end
% 
