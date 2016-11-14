function [sensor] = sensortype(sensorid)
%SENSORTYPE Return sensortype name
%
%  CALL:   sensor = sensortype(sensorid)
%
%  sensorid  = integer defining the sensortype
%  sensor    = string defining the sensortype
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
% validNames = {'n','n_t','n_tt','n_x','n_y','n_xx','n_yy','n_xy','p','u',...
%               'v','w','u_t', 'v_t','w_t','x_p','y_p','z_p','nan'};
% for i=1:19,
%   assert(sensortype(i), validNames{i});
% end
%
% See also sensortypeid, tran

%History
% by pab 2005

  % Note the ordering of validnames can not be changed without changing
  % the order in functions dependent on this function
  validNames = {'n','n_t','n_tt','n_x','n_y','n_xx',...
    'n_yy','n_xy','p','u','v','w','u_t',....
    'v_t','w_t','x_p','y_p','z_p','nan'};
  N = length(validNames)-1;
  if nargin<1||isempty(sensorid)
    sensorid = 1;
  end
  if isnumeric(sensorid)
    sensorid(sensorid<0 | N +1 < sensorid) = N+1;
    sensor = validNames{sensorid};
  else
    error('Input must be an integer!')
  end
  return