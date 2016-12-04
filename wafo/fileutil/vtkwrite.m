function vtkwrite(varargin)
%VTKWRITE Write data to file in VTK-format.
% 
%  CALL: vtkwrite(x,y,z,V,filename,format)
%
%   x,y,z    = vectors or matrices defining the grid for the volume.
%   V        = matrix of volumetric data values
%   filename = name of file  (default tempname) 
%   format   = string defining the VTK format
%              (default 'STRUCTURED_POINTS')
%
% Example: % 
%   
%  [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%  v = x .* exp(-x.^2 - y.^2 - z.^2);
%  vtkwrite(x,y,z,v,'test.vtk')
%  
  
error(nargchk(1,inf,nargin))

validFormats = {'STRUCTURED_POINTS','UNSTRUCTURED_GRID'};
N = nargin;
%Default 
filename  = tempname;
formatStr = 'STRUCTURED_POINTS';

if ischar(varargin{N})
  ind = strmatch(varargin{N},validFormats);
  if ~isempty(ind)
    formatStr = validFormats{ind(1)};
    N = N-1;
  end
end
if ischar(varargin{N})
  filename = varargin{N};
  N = N-1;
else
  %filename = inputname(N);
end  
% make sure the filename is valid 
filename = sprintf('%s.vtk',  strtok(filename,'.'));

if exist(filename,'file') ~= 0, 
  disp(['There allready exist a ', filename ])
  disp(['Copied this file to ', filename  'old'])
  %Make backup copy before overwriting filename
  copyfile(filename,[filename 'old']);
  delete(filename)
end


[fid,msg] =  fopen(filename,'w');

if (fid<0)
  error(msg)
end
commentStr = sprintf('Matlab made file %s',datestr(now));
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'%s\n',commentStr);
fprintf(fid,'ASCII\n\n');

sep=' ';

switch upper(formatStr),
 case 'STRUCTURED_POINTS',
  % reverse the meshgrid order,
  V   = permute(varargin{N},[2 1 3]);
  szV = size(V);
  switch N
   case {1}
    [x0,y0,z0] = deal(0,0,0);
    [dx,dy,dz] = deal(1,1,1);
   case {4}
    [X,Y,Z]    = deal(varargin{1:N-1});
    [x0,y0,z0] = deal(X(1),Y(1),Z(1));
    if ndims(X)>2
      [dx,dy,dz] = deal(diff(X(1,1:2,1)),diff(Y(1:2)),diff(Z(1,1,1:2)));
    else
      [dx,dy,dz] = deal(diff(X(1:2)),diff(Y(1:2)),diff(Z(1:2)));
    end
   otherwise 
    error('1 or 4 numeric inputs expected!')
  end
  
  fprintf(fid,'DATASET %s\n',formatStr);
  fprintf(fid,'DIMENSIONS    %s\n', sprintf('%d  ',szV));
  fprintf(fid,'ORIGIN    %g %g %g\n',x0,y0,z0);
  fprintf(fid,'SPACING    %g %g %g\n\n',dx,dy,dz);
  fprintf(fid,'POINT_DATA   %d\n',prod(szV));
  fprintf(fid,'SCALARS scalars float\n');
  fprintf(fid,'LOOKUP_TABLE default\n\n');
  fprintf(fid,'%g ',V(:));
 case 'UNSTRUCTURED_GRID',
  MX = varargin{N};
  if all(MX(:)==round(MX(:))),
    format = '%d  %d  %d  \n';
    format0 = 'POINTS %d int \n';
  elseif all(MX(:)==round(MX(:)*10)/10),
    format = '%15.2f %15.2f %15.2f  \n';
    format0 = 'POINTS %d float \n';
  else
    format = '%15.4f %15.4f %15.4f  \n';
    format0 = 'POINTS %d float \n';
  end

  fprintf(fid,'%s \n','DATASET UNSTRUCTURED_GRID');
  fprintf(fid,format0,size(MX,1));
  fprintf(fid,format,MX.');
  fprintf(fid,'%s \n',sep);

   
 otherwise
  error('Option is not available.')
end
fclose(fid);