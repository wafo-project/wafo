function out=sphdam(L,n,beta)
%SPHDAM Calculates spherical damage for a 3-D load.
%
%  Calculates the damage on the unit sphere where the damage
%  is defined by
%
%                              b_i
%         D_i(T) =  sum   (x-y)   ,  x>y,
%                  t_j<=T      j
%
%  where  (x,y)_j  is the cycle count counted at time  t_j.
%
%  CALL: D = sphdam(L,n,b);
%
%  where
%
%        D = the damage,
%        L = three column load process,
%        n = the grid size on the unit sphere,
%        b = b_i

data=[ n n beta];

[dimn,dimm]=size(L);
if dimn<dimm, L=L'; end
[dimn,dimm]=size(L);
if dimm~=3
  disp('   Load not tri-axial. Program will terminate.')
  return
end  

disp('   Writing data.')
save sphdam.in data -ascii
save d3load.dat L -ascii

disp('   Starting Fortran executable.')
!sphdam.exe

if nargout==1
  disp('   Loading data.')
  load out.dat
end

