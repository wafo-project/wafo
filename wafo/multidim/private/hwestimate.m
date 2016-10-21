function Hw = HwEstimate(Sf,SfBest,Hw,pos);
%HWESTIMATE  Estimate absolute value of transfer function H(w) from sensor spectra
%
% CALL Hw = hwestimate(Sf,SfBest,Hw,pos)
%
% Internal function to dat2dspec. 
  
%---------------------------------------------------------------------------------
%Estimate the absolute value of the transfer function H(w) from the sensor spectra
%---------------------------------------------------------------------------------
def  = unique(pos(:,4)).';
k0   = find(def==5|def==7|def==11|def==14|def==17);
if any(k0), def(k0) = def(k0)-1;end
k0  = find(def==8);
if any(k0), def(k0) = def(k0)-2;end
def = unique(def);

k = find(SfBest~= 0); % avoid division by zero

for sensorType1 = def
  if (sensorType1==6)  % surface curvatures          : n_xx, n_yy, n_xy
    sensorType2 = 7; sensorType3 = 8;
    ix1 = find(pos(:,4)==sensorType1); Nx1 = length(Nx1);
    ix2 = find(pos(:,4)==sensorType2); Nx2 = length(Nx2);
    ix3 = find(pos(:,4)==sensorType3); Nx3 = length(Nx3);
    Nx  = min([Nx1,Nx2,Nx3]); % need at least one pair of observations
    if (Nx>0)
      Hw0 = sqrt(mean(Sf(ix1(1:Nx),k)+Sf(ix2(1:Nx),k)+2*Sf(ix3(1:Nx),k),1)./SfBest(k));
      %plot(fi(k),Hw0,fi,Hw(ix1(1),:)), legend('Hw0','Hw'), pause
      Hw([ix1;ix2;ix3],k) = Hw0(ones(Nx1+Nx2+Nx3,1),k);
    else 
      warning('Unable to estimate the surface curvature transferfunction!')
    end
  elseif any(sensorType1==[2:3 9 12 15 18])
    ix1 = find(pos(:,4)==sensorType1);
    Nx  = length(ix1);
    if Nx>0
      Hw0       = sqrt(mean(Sf(ix1,k),1)./SfBest(k));
      %plot(fi(k),Hw0,fi,Hw(ix1(1),:)), legend('Hw0','Hw'), pause
      Hw(ix1,k) = Hw0(ones(Nx,1),:);
    else
       warning('Unable to estimate the transferfunction')
    end
  elseif any(sensorType1==[4 10 14 16]) 
    % surface slopes,water particle velocity, water particle acceleration
    % or water particle displacement 
    sensorType2 = sensorType1+1;
    ix1 = find(pos(:,4)==sensorType1); Nx1 = length(ix1);
    ix2 = find(pos(:,4)==sensorType2); Nx2 = length(ix2);
    Nx  = min(Nx1,Nx2); % need at least one pair of observations
    if (Nx>0)
      Hw0 = sqrt(mean(Sf(ix1(1:Nx),k)+Sf(ix2(1:Nx),k),1)./SfBest(k));
      %plot(fi(k),Hw0,fi,Hw(ix1(1),:)), legend('Hw0','Hw'), pause
      
      Hw([ix1;ix2],k) = Hw0(ones(Nx1+Nx2,1),:);
    else %if (any()
      warning('Unable to estimate the transferfunction')
    end
  end
end
return