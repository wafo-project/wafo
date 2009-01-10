function Hs1 = set(Hs,varargin)
%SPECDATA/SET Set SPECDATA object properties
%
%  CALL set(Hs,par1,val1,par2,val2,...)
%
%  par1,par2..= strings identifying the property to alter
%  val1,val2..= corresponding values the parameters are altered to.

error(nargchk(3,inf,nargin));
Hs = parseoptions(Hs,varargin{:});
if ~isa(Hs.wdata,'wdata')
   error('Not allowed to assign a non-WAFO data object to wdata.')
end
if nargout<1
  % Write the newly updated data object back to the calling program
  assignin('caller',inputname(1), Hs);
else
  Hs1 = Hs;
end