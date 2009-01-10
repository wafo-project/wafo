function self = subsasgn(self,index,b,varargin)
%SUBSASGN Define field name assigment for class objects
%
% CALL  self = subsasgn(self,index,b,varargin)
%
% Examples
% f = fdata(1);
% f(2) = fdata(2);
% [f.params] = deal(3); % replace data with 3
% f.params
% f(2).params= 7;
% f.data
% 
% See also set

%History
%By pab 2007

persistent  thisClassTag thisClass
if isempty(thisClassTag)
  thisClass = mfilename('class');
  thisClassTag = sprintf('WAFO:%s:%s',upper(thisClass),upper(mfilename));
end

N = length(index);
subs1 = index(1).subs;

switch index(1).type
case '()'
  if N==1
     if ~isa(b,thisClass)
       error(thisClassTag,'Not allowed to assign a non-WAFO data object to wdata.')
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
%     error(thisClassTag,'???Illegal reference')
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
    error(thisClassTag,'???Illegal reference')
end

