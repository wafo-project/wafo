function SfBest = bfsSpec(Sf,Hw,pos,bfs);
%BFSSPEC Estimate frequency spectrum for the surface elevation from the bfs timeseries
%
%  CALL:  SfBest = bfsSpec(Sf,Hw,pos,bfs);
%
%
def  = unique(pos(bfs,4)).';

k0   = find(def==5|def==7|def==11|def==14|def==17);
if any(k0), def(k0) = def(k0)-1;end

k0  = find(def==8);
if any(k0), def(k0) = def(k0)-2;end

def    = unique(def);
SfBest = Sf;

for sensorType1 = def
  if (sensorType1==6)  % surface curvatures          : n_xx, n_yy, n_xy
    sensorType2 = 7; sensorType3 = 8;
    
    %Find the indices to sensorTypes 1,2 and 3:
    ix1 = find(pos(:,4)==sensorType1); Nx1 = length(Nx1);
    ix2 = find(pos(:,4)==sensorType2); Nx2 = length(Nx2);
    ix3 = find(pos(:,4)==sensorType3); Nx3 = length(Nx3);
    
    Nx  = min([Nx1,Nx2,Nx3]); % need at least one pair of observations
    if (Nx>0)
      k0 = find(Hw(ix1(1),:)==0);
      if any(k0),SfBest([ix1;ix2;ix3],k0)=0 ;end
      
      k   = find(Hw(ix1(1),:)~=0);
      Sf0 = mean(Sf(ix1(1:Nx),k)+Sf(ix2(1:Nx),k)+2*Sf(ix3(1:Nx),k),1)./(Hw(ix1(1),k).^2);
      SfBest([ix1;ix2;ix3],k) = Sf0(ones(Nx1+Nx2+Nx3,1),k);
    else 
      error('Unable to estimate the sea surface spectrum from the surface curvature!')
    end
    
  elseif any(sensorType1==[2:3 9 12 15 18])
    
    % Find indices to sensorType1
    ix1 = find(pos(:,4)==sensorType1);
    Nx  = length(ix1);
    if Nx>0
      k0 = find(Hw(ix1(1),:)==0);
      if any(k0),SfBest(ix1,k0)=0 ;end
      
      k = find(Hw(ix1(1),:)~=0);
      Sf0           = mean(Sf(ix1,k),1)./(Hw(ix1(1),k).^2);
      SfBest(ix1,k) = Sf0(ones(Nx,1),:);
    else
       error('Unable to estimate the sea surface spectrum')
    end
  elseif any(sensorType1==[4 10 14 16]) 
    % surface slopes,water particle velocity, water particle acceleration
    % or water particle displacement 
    sensorType2 = sensorType1+1;
    
    % Find indices to sensorType1 and sensorType2
    ix1 = find(pos(:,4)==sensorType1); Nx1 = length(ix1);
    ix2 = find(pos(:,4)==sensorType2); Nx2 = length(ix2);
    Nx  = min(Nx1,Nx2); % need at least one pair of observations
    if (Nx>0)
      k0 = find(Hw(ix1(1),:)==0);
      if any(k0),SfBest([ix1;ix2],k0)=0 ;end
      
      k   = find(Hw(ix1(1),:)~=0);
      Sf0 = mean(Sf(ix1(1:Nx),k)+Sf(ix2(1:Nx),k),1)./(Hw(ix1(1),k).^2);
      SfBest([ix1;ix2],k) = Sf0(ones(Nx1+Nx2,1),:);
    else %if (any()
      error('Unable to estimate the sea surface spectrum')
    end
  end
end

%Keep only the best frequency spectra
SfBest  = SfBest(bfs,:);
mmx = max(SfBest(:))*1e-5; % Minimum value which is not considered as noise!
k   = find(SfBest < mmx);  % This is most probably noise
if 0, % Geometric mean 
  if any(k)
    SfBest(k) = abs(SfBest(k)+mmx*1e-10);
  end
  if length(bfs) > 1,
    SfBest = exp(mean(log(SfBest))); % geometric mean
  end
else  % Ordinary mean
  if any(k),   % Set the noise to zero. 
    SfBest(k) = 0;
  end
  if length(bfs) > 1,
    SfBest = mean(SfBest);           % ordinary mean
  end
end
return