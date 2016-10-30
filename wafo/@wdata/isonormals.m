function varargout = isonormals(f,p)
%WDATA/ISONORMALS  Isosurface normals.
%
% CALL [fv] = isonormals(f,p)
%
%  f = wdata object
%  p = handle to patch object
%
% Example
% %[x y z v] = flow;
% x = -2:.2:2; y = -2:.25:2; z = -2:.16:2;
% [X,Y,Z] = meshgrid(x,y,z);
% v = X .* exp(-X.^2 - Y.^2 - Z.^2);
% wd = wdata(v,{x,y,z});
% p = patch(isosurface(wd, 0.01));
% isonormals(wd, p)
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
% daspect([1 1 1])
% view(3)
% % camlight; lighting phong
%
% See also isonormals, isosurface

Nf = numel(f);
if Nf>1 
  hold_state = ishold; % remember old hold state
  cfig = gcf;
  
  h = struct('vertices',[],'faces',[]);
  h(Nf).faces  =  [];
  if length(p)==1
    p = p(ones(size(f)));
  end
  for ix=1:Nf,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    h(ix) = isonormals(f(ix),p(ix));
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = isonormals(f.args{:},f.data,p);
end
