function DS = normspfn(DS,thetai)
%NORMSPFN normalizes the spreading function
%         
%  CALL: DS = normspfn(DS,thetai)
%  
% DS     = spreading function size nt x nf
% thetai = angles size nt x 1
% 
% (it also truncate negative values to zero)
%
  
nt  = length(thetai);
ind = find(DS<0 | isnan(DS));
if any(ind),
%  disp('isnan')
%  disp(['Negative directional distribution. Setting negative values to zero. min(DS) = '  num2str(min(DS(ind)))])
  DS(ind) = 0;
end

[ix,iy] = find(DS(:,:) ==inf);
if any(iy)
  for iz = iy(:).' 
    ind0 = (DS(:,iz)<inf); 
    DS(ind0,iz)  = 0;
    DS(~ind0,iz) = 1;
  end
end

%spy(DS),pause

%Normalize so that int D(theta,f) dtheta = 1 for each f 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Sf2      = simpson(thetai,DS);
Sf2 = trapz(thetai,DS);

%plot(Sf2),pause

k        = find(Sf2>sqrt(eps) & Sf2<inf);
if any(k)
  DS(:,k)  = DS(:,k)./Sf2(ones(nt,1),k);
end

return