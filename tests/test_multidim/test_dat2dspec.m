function test_suite=test_dat2dspec()
  initTestSuite;
end
function f = fwaitbar(x, f, msg)
% Trick deliberately overwrite fwaitbar
  if isempty(f),
    f = figure()
  end
end
function test_dat2dspec_()
  S  = jonswap; 
  D  = spreading(linspace(-pi,pi,51),'cos2s'); 
  Sd = mkdspec(S,D,1); 
  Nx = 3; Ny = 2; Nt = 2^14; dx = 10; dy = 10;dt = 0.5; 
  plotflag = 0;
  use_waitbar = 0;
  F = seasim(Sd,Nx,Ny,Nt,dx,dy,dt,1,plotflag, use_waitbar);  
  Z  = permute(F.Z,[3 1 2]); 
  [X,Y] = meshgrid(F.x,F.y); 
  N = Nx*Ny; 
  types = repmat(sensortypeid('n'),N,1); 
  bfs   = ones(N,1); 
  pos   = [X(:),Y(:),zeros(N,1)]; 
  h = inf; 
  Se = dat2dspec([F.t Z(:,:)],[pos types,bfs],h,256,101); % seasim is possibly wrong
end
