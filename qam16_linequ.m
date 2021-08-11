function Lout=qam16_linequ(y,nvar,h,N)

load('signal_constellations');
[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
c=(H*H'+nvar*eye(N))\s;
z=filter(c,1,y);
z=z(1+latency:L+latency);
k=real(c'*s);
Lout=zeros(L*4,1); p=zeros(16,1); l=1;
for n=1:L
    for m=1:16 p(m)=exp((2*real(z(n)*sym_16qam(m)')-k*abs(sym_16qam(m))^2)/(1-k))+1e-10; end
    Lout(l  )=log(p(1)+p(3)+p(5)+p(7)+p(9)+p(11)+p(13)+p(15))-log(p(2)+p( 4)+p( 6)+p( 8)+p(10)+p(12)+p(14)+p(16));
    Lout(l+1)=log(p(1)+p(2)+p(5)+p(6)+p(9)+p(10)+p(13)+p(14))-log(p(3)+p( 4)+p( 7)+p( 8)+p(11)+p(12)+p(15)+p(16));
    Lout(l+2)=log(p(1)+p(2)+p(3)+p(4)+p(9)+p(10)+p(11)+p(12))-log(p(5)+p( 6)+p( 7)+p( 8)+p(13)+p(14)+p(15)+p(16));
    Lout(l+3)=log(p(1)+p(2)+p(3)+p(4)+p(5)+p( 6)+p( 7)+p( 8))-log(p(9)+p(10)+p(11)+p(12)+p(13)+p(14)+p(15)+p(16));
    l=l+4;
end