function [out_llr,v] = polar_decode_L2R(y, u)
N = length(y);
if N == 1        
    out_llr=y;   
    v=u;
else    
    L_w_odd = zeros(1, N/2);
    for index = 1:(N/2)
        L_w_odd(index) = f(y(index), y(index + N/2));
    end    
    [out_llr_1,v_1] = polar_decode_L2R(L_w_odd, u(1:N/2));
    
    L_w_even = zeros(1, N/2);
    for index = 1:(N/2)
        L_w_even(index) = g(y(index), y(index + N/2), v_1(index));
    end
    
    [out_llr_2,v_2] = polar_decode_L2R(L_w_even, u(N/2+1:end));
    
    out_llr=[out_llr_1 out_llr_2; y];  
    
    v = zeros(1, N);
    for index = 1:(N/2)
        v(index) = bitxor(v_1(index), v_2(index));
        v(index + N/2) = v_2(index);
    end    
end
end