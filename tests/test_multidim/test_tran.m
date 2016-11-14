function test_suite=test_tran()
  initTestSuite;
end
function test_tran_()
    N=5000;dt=0.4;f0=0.1;th0=0;h=50; w0 = 2*pi*f0; 
  t = linspace(0,1500,N)'; 
  types = [1 4 5]; 
  bfs   = [1 0 0]; 
  pos   = zeros(3,3); 
  eta0 = exp(-i*w0*t); 
  [Hw Gwt] = tran(w0,th0,pos(1,:),'n',h);  
  eta = real(Hw*Gwt*eta0)+0.001*rand(N,1); 
  [Hw Gwt] = tran(w0,th0,pos(2,:),'n_x',h);  
  n_x = real(Hw*Gwt*eta0)+0.001*rand(N,1); 
  [Hw Gwt] = tran(w0,th0,pos(3,:),'n_y',h);  
  n_y = real(Hw*Gwt*eta0)+0.001*rand(N,1); 
 
  S = dat2dspec([t eta n_x n_y],[pos types',bfs'],h,256,101); 
  plotspec(S) 
 
  close all;
end
