function [NM,H,s,latency]=ut_equinit(h,N)

M=length(h);
NM=N+M-1;
latency=round((NM-1)/2);
H=zeros(N,NM);
for k=1:N for n=k:k+M-1 H(k,n)=h(n-k+1); end; end
u=zeros(NM,1);
u(1+latency)=1;
s=H*u;