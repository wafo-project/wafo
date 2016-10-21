function f_rfc = mc2rfc(f_xy,paramv,paramu)
%MC2RFC  Calculates a rainflow matrix given a Markov chain with kernel f_xy;
%       f_rfc = f_xy + F_mc(f_xy).
%
%  CALL: f_rfc = mc2rfc(f_xy);
%
%  where
%
%        f_rfc = the rainflow matrix,
%        f_xy =  the frequency matrix of Markov chain (X0,X1)
%                but only the triangular part for X1>X0. 
%
%  Further optional input arguments;
%
%  CALL:  f_rfc = mc2rfc(f_xy,paramx,paramy);
%
%       paramx = the parameter matrix defining discretization of x-values,
%       paramy = the parameter matrix defining discretization of y-values,
%      
if nargin<2
paramv=[-1, 1, length(f_xy)];
paramu=paramv;
end

if nargin<3
paramu=paramv;
end
dd=diag(rot90(f_xy));
N=length(f_xy);
Splus=sum(f_xy,2).';
Sminus=fliplr(sum(f_xy));
Max_rfc=zeros(N,1);
Min_rfc=zeros(N,1);
norm=zeros(N,1);
for i=1:N
  Spm=Sminus(i)+Splus(i)-dd(i);
  if Spm>0.
    Max_rfc(i)=(Splus(i)-dd(i))*(Splus(i)-dd(i))/(1-dd(i)/Spm)/Spm;
    Min_rfc(i)=(Sminus(i)-dd(i))*(Sminus(i)-dd(i))/(1-dd(i)/Spm)/Spm;
    norm(i)=Spm;
  end
end

%cross=zeros(N,1);
%for i=2:N
%   cross(N-i+1)=cross(N-i+2)+Sminus(N-i+2)-Splus(N-i+2);
%end

f_rfc=zeros(N,N);
f_rfc(N-1,1)=Max_rfc(N-1);
f_rfc(1,N-1)=Min_rfc(2);

for k=3:N-1
  for i=2:k-1,
    
    %       AAe= f_xy(1:N-k,1:k-i);
    %       SAAe=sum(sum(AAe));
    AA = f_xy(N-k+1:N-k+i,k-i+1:k);
    RAA=f_rfc(N-k+1:N-k+i,k-i+1:k);
    nA=length(AA);
    MA= Splus(N-k+1:N-k+i);
    mA=Sminus(N-k+1:N-k+i);
    normA=norm(N-k+1:N-k+i);
    MA_rfc=Max_rfc(N-k+1:N-k+i);
    %mA_rfc=Min_rfc(k-i+1:k);
    SA=sum(sum(AA));
    SRA=sum(sum(RAA));
    SMA_rfc=sum(MA_rfc);
    SMA=sum(MA);
    DRFC=SA-SMA-SRA+SMA_rfc;
    
    NT=MA_rfc(1)-sum(RAA(1,:));

    %       if k==35
    %          check=[MA_rfc(1) sum(RAA(1,:))]
    %          pause
    %       end

    NT=max(NT,0);

    if NT>1e-6*MA_rfc(1)

      NN=MA-sum(AA,2).';
      e=(fliplr(mA)-sum(AA))';
      e=flipud(e);
      AA=AA+flipud(rot90(AA,-1));
      AA=rot90(AA);
      AA=AA-0.5*diag(diag(AA));
         

      for j=1:nA,
        if normA(j)~=0
          AA(j,:)=AA(j,:)/normA(j);
          e(j)=e(j)/normA(j);
        end
      end
      fx=0.;

      if max(abs(e))>1e-7 && max(abs(NN))>1e-7*MA_rfc(1)
      
 
        I=eye(size(AA));

        if nA==1
          fx=NN/(1-AA)*e;
        else
          fx=NN*((I-AA)\e);
        end
      end
	   
      f_rfc(N-k+1,k-i+1)=DRFC+fx;
    end
  end
  m0=max(0,Min_rfc(N)-sum(f_rfc(N-k+2:N,1)));
  M0=max(0,Max_rfc(N-k+1)-sum(f_rfc(N-k+1,2:k)));
  f_rfc(N-k+1,1)=min(m0,M0);
  %  n_loops_left=N-k+1
end

for k=2:N
  M0=max(0,Max_rfc(1)-sum(f_rfc(1,N-k+2:N)));
  m0=max(0,Min_rfc(k)-sum(f_rfc(2:k,N-k+1)));
  f_rfc(1,N-k+1)=min(m0,M0);
end
f_rfc=f_rfc+rot90(diag(dd),-1);
clf
subplot(1,2,2)
pcolor(levels(paramv),levels(paramu),flipud(f_xy+flipud(rot90(f_xy,-1))))
  axis([paramv(1), paramv(2), paramu(1), paramu(2)])
  title('MC-kernel  f(x,y)')
  ylabel('y'), xlabel('x')                    
axis('square')

subplot(1,2,1)
pcolor(levels(paramv),levels(paramu),flipud(f_rfc))
  axis([paramv(1), paramv(2), paramu(1), paramu(2)])
  title('Rainflow matrix')
  ylabel('max'), xlabel('rfc-min')                    
axis('square')








