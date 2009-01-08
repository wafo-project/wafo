function out = set(self,varargin)
%SPECDATA/SET Set SPECDATA object properties
%
%  CALL set(Hs,par1,val1,par2,val2,...)
%
%  par1,par2..= strings identifying the property to alter
%  val1,val2..= corresponding values the parameters are altered to.


ni = nargin;
no = nargout;

if ~isa(self,'specdata'),
   % Call built-in SET. Handles calls like set(gcf,'user',ss)
   builtin('set',self,varargin{:});
   return
end
if ni==1,
 display(self)
else
  self =  parseoptions(self,varargin{:});
  if ~isa(self.wdata,'wdata')
    error('Not allowed to assign a non-WAFO data object to wdata.')
  end

  if no<1
    name = inputname(1);
    if ~isempty(name)
      % Write the newly updated data object back to the calling program
      assignin('caller',inputname(1),self);
    else
      warning('WAFO:SPECDATA:SET','Unable to update SPECDATA object, because the input argument is a result of a calculation or an expression!')
    end
  else
    out =  self;
  end
end