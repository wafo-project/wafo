function NT = cmat2nt(F,def)
% CMAT2NT Calculates a counting distribution from a cycle matrix.
%
% CALL:  NT = cmat2nt(F,def);
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
%   F = round(triu(rand(4),1)*10)
%   NT = cmat2nt(F)
%
% See also  nt2cmat

% Tested on Matlab 6.0
%
% History:
% Created by PJ (Pär Johannesson) 19-Nov-1999
% Earlier version named 'nt2fr' in WAT

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni<2
  def = 1;
end

n=length(F); % Number of discrete levels

if def == 1
  
  NT = zeros(n);
  NT(2:n,1:n-1) = fliplr(cumsum(fliplr(cumsum(F(1:end-1,2:end))),2));
  
elseif def == 11 % same as def=1 but using for-loop
  
  NT = zeros(n);
  for i = 2:n
    for j= 1:n-1
      NT(i,j) = sum(sum(F(1:i-1,j+1:n)));
    end
  end
  
elseif def == 0
  
  disp(['def = ' num2str(def) ' not yet implemented'])
  
elseif def == -1
  
  disp(['def = ' num2str(def) ' not yet implemented'])
  
else
  
  disp(['def = ' num2str(def) ': not a valid value of def'])
  
end


  
