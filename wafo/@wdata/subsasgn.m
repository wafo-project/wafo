function self = subsasgn(self,index,b,varargin)
%WDATA/SUBSASGN Define field name assigment for WDATA objects
%
% CALL  self = subsasgn(self,index,b,varargin)
%
% Examples
% f = wdata(1);
% f(2) = wdata(2);
% [f.data] = deal(3); % replace data with 3
% f.data
% f(2).data= 7;
% f.data
% 
% See also wdata/set

%History
%By pab 2007

N = length(index);
subs1 = index(1).subs;

switch index(1).type
case '()'
  if N==1
     if ~isa(b,'wdata')
       error('WAFO:WDATA:SUBSASGN','Not allowed to assign a non-WAFO data object to wdata.')
     end
    self(subs1{:}) = b;
  else
     self(subs1{:}) = subsasgn(self(subs1{:}),index(2:N),b,varargin{:});
%     if strcmpi(index(2).type,'.')
%     if N==2,
%       [self(subs1{:}).(index(2).subs)] = deal(b,varargin{:});
%     else
%       pself = self(subs1{:});
%        par = cell(1,numel(pself));
%       [par{:}] = deal(b,varargin{:});
%       for ix = 1:numel(pself)
%         pself(ix).(subs1) = subsasgn( self(ix).(subs1),index(3:N),par{ix});
%       end
%       self(subs1{:}) = pself;
%     end
%   else
%     error('WAFO:WDATA:SUBSASGN','???Illegal reference')
  end
  
  case '.'
    if N==1
      [self.(subs1)] = deal(b,varargin{:});
    else
      par = cell(1,numel(self));
      [par{:}] = deal(b,varargin{:});
      for ix = 1:numel(self)
        self(ix).(subs1) = subsasgn( self(ix).(subs1),index(2:N),par{ix});
      end
    end
  otherwise
    error('WAFO:WDATA:SUBSASGN','???Illegal reference')
end

