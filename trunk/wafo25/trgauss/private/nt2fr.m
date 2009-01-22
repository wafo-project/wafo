function F=nt2fr(M,def)
%NT2FR  Calculates the frequency matrix given the counting distribution matrix.
%
%  CALL: fr = nt2fr(NT);
%        fr = nt2fr(NT,def);
%
%  where
%
%        fr  = the frequency matrix,
%        NT  = the counting distribution matrix,
%        def = 0,1
%        k   = number of diagonals from main diagonal to be
%              set to zero (optional input argument).
%
%  If def=0 function computes the inverse to 
%
%     N_T(u,v) = #{ (M_i,m_i); M_i>u, m_i<v }
%
%  and if def=1 the inverse to
%
%     N_T(u,v) = #{ (M_i,m_i); M_i>=u, m_i=<v }
%
%  where  (M_i,m_i)  are  cycles and v,u are in the discretization grid.
%


n=length(M); P=zeros(size(M));

if nargin<2
def=0;
end

if def==0
P(1:n-1,1:n-1)=M(1:n-1,1:n-1)+M(2:n,2:n)-M(1:n-1,2:n)-M(2:n,1:n-1);
else
P(2:n,2:n)=M(1:n-1,1:n-1)+M(2:n,2:n)-M(1:n-1,2:n)-M(2:n,1:n-1);
end

%if nargin==2
%  k_cut=k;
%else
%  k_cut=2;
%end  

F=P;

%F=fliplr(triu(fliplr(P),k_cut));


