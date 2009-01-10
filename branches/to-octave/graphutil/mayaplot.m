function mayaplot(varargin)
%MAYAPLOT Visualize the data using MayaVi
%
% CALL:  mayaplot(X,Y,Z,V,filename)
%        mayaplot(V,filename)  
%
%        X,Y,Z = matrices defining the grid for the volume.
%        V     = Volumetric data to vizualize 
%     filename = string giving the name of the file, possibly with the
%                complete path. (default tempname)  
%
% MAYAPLOT Slice and isosurface volume exploration GUI
%
% Installation note:  
% Download and install MayaVi1 from:
%
% <http://mayavi.sourceforge.net/>
%
% After installing mayavi, add the directory where mayvi.exe
% is located into the path variable (environment variable).
% That is all. Now you are ready to visualize your data.
% 
% Example:
%  [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%  v = x .* exp(-x.^2 - y.^2 - z.^2);
%  mayaplot(x,y,z,v)
%
% See also: vtkwrite  

  
   
%Tested on: Matlab 6.1  Matlab 5.3
% revised pab May 2005
% -fixed a bug: filename now works correctly
% revised pab June 2004
% -added filename  
%By Per Andreas Brodtkorb 02.04.2003

% Installation note:  
% After installing mayavi, add the directory where mayvi.exe
% is located into the path variable. That is all. Now you are ready to
% visualize your data.
  
  
error(nargchk(1,5,nargin))
np = nargin;
if ischar(varargin{np})
  filename = strtok(varargin{np},'.'); % remove extension if any
  np = np-1;
else
  filename = tempname;
  
end

vtkfilename = sprintf('%s.vtk',filename);
vtkwrite(varargin{1:np},vtkfilename);

exeStr = sprintf('mayavi -d %s -m Axes -m GridPlane -M new -f Threshold -m IsoSurface &', vtkfilename);
s = dos(exeStr);
if s==0
  error('MayaVi could not be found. Check that it is installed?')
end
return