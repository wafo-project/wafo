function r=fr2res(f)
%FR2RES Generates a stationary residual from the frequency matrix.

%  Copyright 1993, Mats Frendahl & Igor Rychlik,
%  Dept. of Math. Stat., University of Lund.

% Changed by PJ (Pär Johannesson) 16-Feb-2004

N=length(f); 

% Get index for upper left element > 0 in frequency matrix.

[index_i,index_j]=getulcc(f);

if N+1-index_j-index_i < 4 
    r=[N+1-index_j index_i];
    return
end

left_done=0; right_done=0; success=0; direction=-1;

while ~success
    
  % Let the residual be (min,max) with levels found above.
  
  r=[N+1-index_j index_i];
    
  comb=fr2comb(f,r);
    
  X=r(1); again=1; number_of_times=0;
    
  while again
        
    h=0;
        
    if r(1)>r(2)
      % disp('upslope')
            
      if r(1)>(r(2)+1)
        m=r(1); M=r(2); j=N-m+1; i=m-1:-1:M+1;
        combinations=comb(i,j);
        freq=f(i,j);
        prob=comb2pro(combinations,freq);
        h=min(find(rand<cumsum(prob)))-1;
      end
            
      if h>0
        C=m-h;
        i=C:m;
        if f(C,j)>0
          r=[m-1 C r];
        else
          warning('Warning 1')
          break
          %r=[r(1)-1 r(2:length(r))];
        end
      else
        r(1)=r(1)-1;
      end
      
    else
      
      % disp('downslope')
      if r(1)<r(2)-1
        m=r(2); M=r(1); i=M; j=M+1:m-1; j=N-j+1;
        combinations=comb(i,j)';
        freq=f(i,j)';
        prob=comb2pro(combinations,freq);
        h=min(find(rand<cumsum(prob)))-1;
      end
      
      if h>0
        j=M:M+h; j=N-j+1;
        c0=M+h;
        c=N+1-c0;
        if f(i,c)>0
          r=[M+1 c0 r];
        else
          warning('Warning 2')
          break
          %r=[r(1)+1 r(2:length(r))];
        end
      else
        r(1)=r(1)+1;
      end
    end
    
    if r(1)==r(2), r=r(2:length(r)); end
    
    %if (length(r)>=2), again=1; else again=0; end
    
    X=[X r(1)];
    
    number_of_times=number_of_times+1;
    
    if length(r)<2
      r=[N+1-index_j index_i];
      number_of_times=0;
    end
    
    again=number_of_times<100;
    
  end
  
  if ~isempty(r)
    if direction==-1
      r1=r;
      left_done=1;
      direction=1;
    elseif direction==1
      r2=r;
      right_done=1;
    end
  end
  
  success=left_done&right_done;
  
end

r=[r1(1:length(r1)) N+1-index_j  fliplr(r2)];

%subplot(2,2,1)
%plot(r1)
%subplot(2,2,2)
%plot(fliplr(r2))
%subplot(2,2,3)
%plot(r)

