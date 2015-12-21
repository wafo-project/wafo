function F = nt2cmat(NT,def)
% NT2CMAT  Calculates a cycle matrix from a counting distribution.
%
% CALL:  F = nt2cmat(NT,def);
%
%   NT       = Counting distribution. [nxn]
%
%   F        = Cycle matrix. [nxn]
%   def      =  1: causes peaks to be projected upwards and troughs 
%                  downwards to the closest discrete level (default).
%            =  0: causes peaks and troughs to be projected to 
%                  the closest discrete level.
%            = -1: causes peaks to be projected downwards and the 
%                  troughs upwards to the closest discrete level.
%
% Example: 
%   F0 = round(triu(rand(4),1)*10)
%   NT = cmat2nt(F0)
%   F = nt2cmat(NT)
%
% See also  cmat2nt

% Tested on Matlab 6.0
%
% History:
% Revised by PJ 18-May-2000
%   Updated help text.
% Created by PJ (Pär Johannesson) 23-Nov-1999
% Earlier version named 'nt2fr' in WAT
  
% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni<2
  def = 1;
end

n=length(NT); % Number of discrete levels

if def == 1
  
  F = zeros(n);
  I=1:n-1;
  J=2:n;
  F(I,J) = NT(I+1,J-1)-NT(I,J-1)-NT(I+1,J)+NT(I,J);
  
elseif def == 11 % same as def=1 but using for-loop
  
  F = zeros(n);
  for i = 1:n-1
    for j= 2:n
      F(i,j) = NT(i+1,j-1)-NT(i,j-1)-NT(i+1,j)+NT(i,j);
    end
  end
  
elseif def == 0
  
  disp(['def = ' num2str(def) ' not yet implemented'])
  
elseif def == -1
  
  disp(['def = ' num2str(def) ' not yet implemented'])
  
else
  
  disp(['def = ' num2str(def) ': not a valid value of def'])
  
end
