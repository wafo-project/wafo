function  self = settype(self)
%SETTYPE Determine the data type data_1d, data_2d, data_3d or data_nd
%
% CALL: self = settype(self)

Nf = numel(self);
if Nf>1
  for ix = 1:Nf
    self(ix) = settype(self(ix));
  end  
else
  if iscell(self.args) % Multidimensional data
    switch length(self.args);
      case {0,1},
        txt = sprintf('%s\n','Unable to determine type, because length(self.args)<2.',...
          'If the data is 1D, then self.args should be a vector!',...
          'If the data is 2D, then length(self.args) should be 2.',...
          'If the data is 3D, then length(self.args) should be 3.',...
          'Unless you fix this, the plot methods will not work!');
        Warning('WAFO:WDATA:SETTYPE',txt)
      case 2,
        self.type = data_2d;
      case 3,
        self.type = data_3d;
      otherwise
        self.type = data_nd;
    end
  else % One dimensional data
    self.type = data_1d;
  end
end % 

if nargout<1
  % Write the newly updated data object back to the calling program
  assignin('caller',inputname(1), self);
end
