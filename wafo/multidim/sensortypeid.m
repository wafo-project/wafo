function sensorid = sensortypeid(sensortype)
%SENSORTYPEID Return ID for sensortype name
%
%  CALL:   sensorid = sensortypeid(sensortype)
%
%  sensorid   = integer defining the sensortype
%  sensortype = string defining the sensortype
%               valid options are:
%   'n'    : Surface elevation              (n=Eta)  (default)
%   'n_t'  : Vertical surface velocity
%   'n_tt' : Vertical surface acceleration
%   'n_x'  : Surface slope in x-direction 
%   'n_y'  : Surface slope in y-direction
%   'n_xx' : Surface curvature in x-direction
%   'n_yy' : Surface curvature in y-direction
%   'n_xy' : Surface curvature in xy-direction
%   'P'    : Pressure fluctuation about static MWL pressure 
%   'U'    : Water particle velocity in x-direction
%   'V'    : Water particle velocity in y-direction
%   'W'    : Water particle velocity in z-direction
%   'U_t'  : Water particle acceleration in x-direction
%   'V_t'  : Water particle acceleration in y-direction
%   'W_t'  : Water particle acceleration in z-direction
%   'X_p'  : Water particle displacement in x-direction from its mean position
%   'Y_p'  : Water particle displacement in y-direction from its mean position
%   'Z_p'  : Water particle displacement in z-direction from its mean position
%
% Example:
% assert(sensortypeid(strvcat('W','v', 'n')), [12,11,1]')
% assert(isnan(sensortypeid('rubbish')))
%
% See also sensortype, id

%History
% by pab 2005

  % Note the ordering of validnames can not be changed without changing
  % the order in functions dependent on this function
  validNames = strvcat('n','n_t','n_tt','n_x','n_y','n_xx',...
                       'n_yy','n_xy','p','u','v','w','u_t',....
                       'v_t','w_t','x_p','y_p','z_p');
  if nargin<1||isempty(sensortype)
    sensortype = 'n';
  end
  if ischar(sensortype)
    %[C,ia,ib] = intersect(validNames,lower(sensortype),'rows');
    N1 = size(sensortype,1);
    sensorid = repmat(nan,N1,1);
    
    for ix = 1:N1
      tmp =    strmatch(lower(deblank(sensortype(ix,:))),validNames,'exact');
      if ~isempty(tmp)
        sensorid(ix) = tmp;
      end
    end
  elseif isnumeric(sensortype)
    N = size(validNames,1);
    sensorid = sensortype;
    sensorid(sensorid>N) = nan;
  else 
    error('Input must be character arrays!')
  end
  return