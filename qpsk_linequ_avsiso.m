function Lout=qpsk_linequ_avsiso(y,nvar,h,N,Lin)

[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
th=tanh(Lin*0.5);
xm=(th(1:2:2*L)+th(2:2:2*L)*i)/sqrt(2);
av=mean(1-xm.*conj(xm));
xm=[xm;zeros(N,1)];
y=y-filter(h,1,xm);
c=(H*H'*av+(1-av)*s*s'+nvar*eye(N))\s;
k=real(c'*s);
p=real(H'*c);
q=real(nvar*c'*c);
z=filter(c,1,y);
z=z(1+latency:L+latency);
z=z+k*xm(1:L);
zv=zeros(L,1);
for n=1+latency:L+latency
    if (n<NM) 
        x=flipud([zeros(NM-n,1);xm(1:n)]);
        %Downward circular shift for a matrix
        %Transpose and order reverse for a column vector
    else
        x=flipud(xm(n-NM+1:n));
    end
    x(1+latency)=1; 
    zv(n-latency)=p'*diag(1-x.*conj(x))*p+q;
end
Lout=zeros(2*L,1);
Lout(1:2:2*L)=sqrt(8)*k*real(z)./zv;
Lout(2:2:2*L)=sqrt(8)*k*imag(z)./zv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: LLRs calculated with constant variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lout(1:2:2*L)=sqrt(8)*k*real(z)/(q+av*(p'*p-k'*k));
% Lout(2:2:2*L)=sqrt(8)*k*imag(z))/(q+av*(p'*p-k'*k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






