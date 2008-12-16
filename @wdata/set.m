function out = set(self,varargin)
%WDATA/SET Set WDATA object properties
%
%  CALL Hs = set(Hs,par1,val1,par2,val2,...)
%
%  par1,par2..= strings identifying the property to alter
%  val1,val2..= corresponding values the parameters are altered to.
%
% Examples
% d = wdata;
% set(d)                 % displays all properties of d and their values.  
%
% set(d,'data',1)        % equals: d.data = 1
% d.data
%
% d1 = set(d,'data',2)   % equals: d1 = d; d1.data = 2
% d2 = set(d,'data',3,'args',1) % equals d2 = d; d2.data=3;d2.args = 1;
%
% d3 = [d1,d2];
% d3(1).set('data',5) % equals d3(1).data=5;
% d3.data
%
% set(d3(2),'data',5) % is illegal call. d3 is not changed.
% d3.data
% d3(2) = set(d3(2),'data',5) % But this is OK
% d3.data
% 
% See also wdata/get


% History
% revised pab april 2007
% -added warning if input object is a result of a computation
% By pab 2007

ni = nargin;
no = nargout;

% if ~isa(self,'wdata'),
%    % Call built-in SET. Handles calls like set(gcf,'user',ss)
%    builtin('set',self,varargin{:});
%    return
% end


% Handle read-only cases
if ni==1,
 display(self)
elseif no<1
  % Write the newly updated data object back to the calling program
  name = inputname(1);
  if ~isempty(name)
    assignin('caller',name, parseoptions(self,varargin{:}));
  else
    warning('WAFO:WDATA:SET','Unable to update WDATA object, because the input argument is a result of a calculation or an expression!')
  end
else
  out = parseoptions(self,varargin{:});
end

