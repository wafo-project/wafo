function f_rfc = mctp2rfc(f_mM,f_Mm)
%MCTP2RFC  Rainflow matrix given a Markov matrix of a Markov chain of turning points
%
%  computes f_rfc = f_mM + F_mct(f_mM).
%
%  CALL: f_rfc = mctp2rfc(f_mM);
%
%  where
%
%        f_rfc = the rainflow matrix,
%        f_mM =  the min2max Markov matrix, 
%
%  Further optional input arguments;
%
%  CALL:  f_rfc = mctp2rfc(f_mM,f_Mm,paramm,paramM);
%
%        f_Mm  = the max2min Markov matrix, 
%       paramm = the parameter matrix defining discretization of minimas,
%       paramM = the parameter matrix defining discretization of maximas,
% 
% Example
%  fmm = [0.0183    0.0160    0.0002    0.0000         0;...
%    0.0178    0.5405    0.0952         0         0;...
%    0.0002    0.0813         0         0         0;...
%    0.0000         0         0         0         0;...
%         0         0         0         0         0]
%  mctp2rfc(fmm)

%%
% TODO %  Check this: paramm and paramM are never used?????

if nargin==1
  f_Mm=f_mM;
end
      
% if nargin<3
%   paramm=[-1, 1 ,length(f_mM)];
%   paramM=paramm;
% end
% 
% if nargin<4
%   paramM=paramm;
% end


N=length(f_mM);
Max=sum(f_mM,2).';
Min=sum(f_mM);
f_rfc=zeros(N,N);
f_rfc(N-1,1)=Max(N-1);
f_rfc(1,N-1)=Min(N-1);
for k=3:N-1
  for i=2:k-1,
    AA=f_mM(N-k+1:N-k+i-1,k-i+1:k-1);
    AA1=f_Mm(N-k+1:N-k+i-1,k-i+1:k-1);
    RAA=f_rfc(N-k+1:N-k+i-1,k-i+1:k-1);
    nA=length(AA);
    MA=Max(N-k+1:N-k+i-1);
    mA=Min(k-i+1:k-1);
    SA=sum(sum(AA));
    SRA=sum(sum(RAA));
    %SMA=sum(MA);
    %SmA=sum(mA);
    DRFC=SA-SRA;
    NT=min(mA(1)-sum(RAA(:,1)),MA(1)-sum(RAA(1,:)));
    NT=max(NT,0);
    
    if NT>1e-6*max(MA(1),mA(1))

      NN=MA-sum(AA,2).';
      e=(mA-sum(AA))';
      e=flipud(e);
      PmM=rot90(AA);
      for j=1:nA,
        norm=mA(nA-j+1);
        if norm~=0
          PmM(j,:)=PmM(j,:)/norm;
          e(j)=e(j)/norm;
        end
      end
      fx=0.;

      if max(abs(e))>1e-6 && max(abs(NN))>1e-6*max(MA(1),mA(1))
      
        PMm=AA1;
        for j=1:nA,
          norm=MA(j);
          if norm~=0
            PMm(j,:)=PMm(j,:)/norm;
          end
        end
        PMm=fliplr(PMm);
 
        A=PMm; B=PmM;
        I=eye(size(A));

        if nA==1
          fx=NN*(A/(1-B*A)*e);
        else
         fx=NN*(A*((I-B*A)\e));
        end
      end
	   
      f_rfc(N-k+1,k-i+1)=fx+DRFC;

      %  check2=[ DRFC  fx]
      % pause
    else
      f_rfc(N-k+1,k-i+1)=0.;
    end
  end
  m0=max(0,Min(1)-sum(f_rfc(N-k+2:N,1)));
  M0=max(0,Max(N-k+1)-sum(f_rfc(N-k+1,2:k)));
  f_rfc(N-k+1,1)=min(m0,M0);
  %  n_loops_left=N-k+1
end

for k=2:N
  M0=max(0,Max(1)-sum(f_rfc(1,N-k+2:N)));
  m0=max(0,Min(N-k+1)-sum(f_rfc(2:k,N-k+1)));
  f_rfc(1,N-k+1)=min(m0,M0);
end

%clf
%subplot(1,2,2)
%pcolor(levels(paramm),levels(paramM),flipud(f_mM))
%  title('Markov matrix')
%  ylabel('max'), xlabel('min')                    
%axis([paramm(1) paramm(2) paramM(1) paramM(2)])
%axis('square')

%subplot(1,2,1)
%pcolor(levels(paramm),levels(paramM),flipud(f_rfc))
%  title('Rainflow matrix')
%  ylabel('max'), xlabel('rfc-min')                    
%axis([paramm(1) paramm(2) paramM(1) paramM(2)])
%axis('square')











