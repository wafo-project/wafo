function varargout = isosurface(f,varargin)
%WDATA/ISOSURFACE  Isosurface extractor.
%
% CALL [fv] = isosurface(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see isosurface for details.
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
% % alpha(0.5) % set transparency
%
% See also isosurface

  Nf = numel(f);
  if Nf>1 
    hold_state = ishold; % remember old hold state
    cfig = gcf;
    
    h = struct('vertices',[],'faces',[]);
    h(Nf).faces  =  [];
    for ix=1:Nf,
      if hold_state
        newplot
      else
        figure(cfig-1+ix)
      end
      h(ix) = isosurface(f(ix),varargin{:});
    end  
    if nargout>0
      varargout{1} = h;
    end
  else
     args = _mesh_args(f);
     [varargout{1:nargout}] = isosurface(args{:}, f.data, varargin{:});
     labelfig(f);
  end
end

