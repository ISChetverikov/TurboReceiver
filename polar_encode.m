function x = polar_encode(u,N)
    n = log2(N);
    m = 1; %number of bits combined
    for d = n-1:-1:0
        for i= 1:2*m:N
            a = u(i:i+m-1); %first part
            b = u(i+m:i+2*m-1); %second part
            u(i:i+2*m-1) = [mod(a+b,2) b]; %combining
        end
        m = m * 2;
    end
    x = u;
end