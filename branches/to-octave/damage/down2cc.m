function [cc,delta] = down2cc(mu,n)
%DOWN2CC Calculates the cycle count which has the highest damage 
%  given the downcrossing intensity. 
%
% CALL:  [cc,delta] = down2cc(cross,n);
%
%       cc    = the cycle count,
%       delta = the slice distance,
%       cross = a two column matrix with the levels  u  in the first
%               column and the number of downcrossings/downcrossing 
%               intensity in the second,
%       n     = the number of slice levels between 0 and maximum of
%               the number of downcrossings/downcrossing intensity.



%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

% History
% Revised pab 2007
% -translated fortran code into matlab which is faster because of
% the nececity of writing to the harddisk.
% -moved code calling fortran version into a subfunction
% By Mats Frendahl 1993

try
  [cc,delta] = down2cc_m(mu,n);
catch
  [cc,delta] = down2cc_f77(mu,n);  
end




function [cc,delta] = down2cc_f77(mu,n)
%DOWN2CC_F77 fortran version


filename1='slice.dat';
if exist(filename1,'file')
  delete(filename1)
end
fid = fopen(filename1,'wt');
fprintf(fid,'%4.0f \n',n);
fclose(fid);

filename2='crossint.dat';
if exist(filename2,'file')
  delete(filename2)
end
fid = fopen(filename2,'wt');
fprintf(fid,'%10.5f %10.5f\n',mu');
fclose(fid);


disp('   Starting Fortran executable.')
dos([ wafoexepath 'down2cc.exe']);

disp('   Loading data.')

out = load('out.dat');
cc=out(:,2:3);

delta=max(mu(:,2))/n;

delete(filename1,filename2,'out.dat')


function  [cc,delta] = down2cc_m(mu,num_slices)
%DOWN2CC_M matlab version

%     Initiate, max_mu, the maximum value
%     of the crossing intensity.
max_mu = max(mu(:,2));
%     The slice distance is calculated.
delta  = max_mu/num_slices;
      

u = mu(:,1);

cc = zeros(10*num_slices,2);      
%     Initiate, nc, the number of cycles. 
nc=0;

%     Fix a slice level.
      
for ii=1:num_slices

  %     Calculate the height of the slice level. OBS: We lower the
  %     slice levels by delta/2 in order to a) really get a meaning-
  %     full cycle if the number of slice levels = 1 and b) to get
  %     closer to the bottom where the "tail-cycles" can be hard to get.
         
  h=(ii-0.5)*delta;


  %     Run trough all  u  and find where the crossing intensity
  %     crosses the slice level.
  
  % Note this gives slightly different results
  % because inflection points are included in the fortran version of 
  % down2cc but not in findcross used here
  ind = findcross(mu(:,2),h);   
         
  n = length(ind);
  %     If the number of crossings are greater than 0, extract
  %     the cycles by combining a cycle using the two closest
  %     u's that cross the slice level. We take the mean of the
  %     two u's near the crossing. We also save the slice level.
         
  
  for k=1:2:n
    nc=nc+1;
    %out(nc,1)=h
    cc(nc,1)=(u(ind(k+1))+u(ind(k+1)+1))/2;
    cc(nc,2)=(u(ind(k))+u(ind(k)+1))/2;
  end
end

cc(nc+1:end,:) = [];
