function [b,varargout] = subsref(self,index)
%WDATA/SUBSREF Define field name indexing for objects or methods to apply
%
% 
% Examples
%  x = linspace(0,5); 
%  d = wdata(pdfray(x,1),x);
%  d.plot               % plot contents of d
%  d.argMax
%
%  d2 = d; d2(2) = wdata(pdfray(x,0.5),x);
%  am = {d2.argMax}            % cellarray of argMax values
%  am1 = {d2(1:2).argMax}      % NB!: Matlabbug returns only the first value
%  [am2{1:2}] = d2(1:2).argMax % Use this call instead!
%
%  d2.plot;           % plot contents of d2 in separate figures
%  d2(1).plot         % plot only 1'st object only.
%  hold on
%  d2.plot;           % plot both objects in the same figure
%  
%
%
% See also wdata/get

% TODO % This is slow. Optimize this for speed.
% TODO % am1 = {d2(1:2).argMax} is a Matlabbug returns only the first value

%History
%By pab 2007


persistent legalMethodName legalFieldName outPutMethods
if isempty(legalMethodName) || isempty(legalFieldName) || isempty(outPutMethods)
  legalMethodName = methods(self,'-full');
  legalFieldName = fieldnames(self);
  outPutMethods = {'char','fieldnames'};
end

N = length(index);

varargout = cell(1,nargout-1);

%if 1
switch index(1).type
  case '()'
    % Working only on a part of self
    pself  = self(index(1).subs{:}); % copy part of self
    if N == 1
      b =  pself;
    elseif strcmpi(index(2).type,'.')
      switch index(2).subs
        case legalMethodName
          N3 = min(N,3);
          if nargout<1 && ~any(strcmpi(index(2).subs,outPutMethods))
            feval(index(2).subs,pself,index(3:N3).subs{:});
            if strncmpi(index(2).subs,'set',3)
              self(index(1).subs{:}) = pself;
              assignin('caller',inputname(1), self);
            end
          else
            nout = min(nargout(index(2).subs),nargout);
            [b,varargout{1:nout-1}] = feval(index(2).subs,pself,index(3:N3).subs{:});
          end
        otherwise
          %nout = max(min(numel(pself),nargout),1);
          nout = numel(pself);
          varargout = cell(1,nout-1);
          [b, varargout{1:nout-1}] = subsref(pself,index(2:end));
      end % switch
    else
      error('WAFO:WDATA:SUBSREF','???? Reference to field or method expected.')
    end % if
    
    
  case '{}'
    error('WAFO:WDATA','??? Cell contents reference from a non-cell array object.')
    
    
  case '.'
    Nself = numel(self);
    switch index(1).subs
      case legalFieldName
        nout = max(min(Nself,nargout),1);
        if N == 1
          [b,varargout{1:nout-1}] =  deal(self(1:nout).(index.subs));
        else          
          %[b,varargout{1:nargout-1}] =  deal(subsref(self.(index(1).subs),index(2:end)));
           b = subsref(self(1).(index(1).subs),index(2:end));          
            for ix=2:nout
              varargout{ix-1} =  subsref(self(ix).(index(1).subs),index(2:end));
            end
            if nout<2 &&  nargout>1
              [varargout{:}] = deal(b);
            end
        end
      case legalMethodName
        N2 = min(N,2);
        if nargout<1 && ~any(strcmpi(index(1).subs,outPutMethods))
          feval(index(1).subs,self,index(2:N2).subs{:});
          if strncmpi(index(1).subs,'set',3)
            assignin('caller',inputname(1), self);
          end
        else
          nout = min(nargout(index(1).subs),nargout);
          [b,varargout{1:nout-1}] = feval(index(1).subs,self,index(2:N2).subs{:});
         
        end
      otherwise
        error('WAFO:WDATA','???? Reference to non-existent field or method %s',index(1).subs)
    end % switch
end % switch
  

return
%end

% Old call Kept just in case
% b = self;
% for ix = 1:N
%   subs = index(ix).subs;
%   switch index(ix).type
%     case '()'
%      b =  b(subs{:});
%     case '{}'
%       if ix == N
%         [b,varargout{1:nargout-1}] = deal(b{subs{:}});
%         return
%       else
%         b =  b{subs{:}};
%       end
%     case '.'
%        if ix == N
%           [b,varargout{1:nargout-1}] = deal(b.(subs));
%           return
%        elseif (ix == 1) || (ix==2)
%           if any(strcmpi(subs,legalFieldName))
%             b = b.(subs);
%           elseif (N==2 || N==3) && strcmpi(index(N).type,'()')
%             
%             iz = find(strcmpi(subs,legalMethodName));
%             if any(iz)
%               if nargout<1 && ~any(strcmpi(legalMethodName{iz},outPutMethods))
%                 feval(legalMethodName{iz},b,index(ix+1).subs{:});
%                 if isa(b,'wdata')
%                   assignin('caller',inputname(1), b);
%                 end
%                 clear b
%               else
%                 [b,varargout{1:nargout-1}]  = feval(legalMethodName{iz},b,index(ix+1).subs{:});
%               end
%             else
%               error('???? Reference to non-existent field or method %s',subs)
%             end
%             return
%           else
%             error('???? Reference to non-existent field %s',subs)
%           end
%          
%        else
%           b = b.(subs);
%        end
%   end
% end
% return
% 
% 
