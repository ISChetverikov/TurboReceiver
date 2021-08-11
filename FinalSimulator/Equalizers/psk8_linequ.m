function Lout=psk8_linequ(y,nvar,h,N)

load('signal_constellations');
[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
c=(H*H'+nvar*eye(N))\s;
z=filter(c,1,y);
z=z(1+latency:L+latency);
k=real(c'*s);
Lout=zeros(L*3,1); p=zeros(8,1); l=1;
for n=1:L
    for m=1:8 p(m)=exp(2*real(z(n)*sym_8psk(m)')/(1-k))+1e-10; end
    Lout(l  )=log(p(1)+p(3)+p(5)+p(7))-log(p(2)+p(4)+p(6)+p(8));
    Lout(l+1)=log(p(1)+p(2)+p(5)+p(6))-log(p(3)+p(4)+p(7)+p(8));
    Lout(l+2)=log(p(1)+p(2)+p(3)+p(4))-log(p(5)+p(6)+p(7)+p(8));
    l=l+3;
end