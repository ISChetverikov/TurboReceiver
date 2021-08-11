function L = MseEqualize(y, H, sigma2, xm)
    pskDemodulator = comm.PSKDemodulator( ...
                'ModulationOrder',  2^2, ...
                'PhaseOffset',      pi/4, ...
                'BitOutput',        true, ...
                'Variance',         sigma2, ...
                'DecisionMethod', 'Log-likelihood ratio');
            
    N = length(y); 
    L = zeros(2 * N, 1);
    c = (H'*H + sigma2) \ 1;
    k = real(c * 1);
    y1 = y;
    if (~isempty(xm))
        y1 = y + xm;
    end
    z = c * H' * y1;
    % L(1:2:2*N) = sqrt(8) * real(z) / (1-k);
    % L(2:2:2*N) = sqrt(8) * imag(z) / (1-k);
    L = pskDemodulator(z');
end