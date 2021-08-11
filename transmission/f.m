function z = f(x, y)
    z = sign(x)*sign(y)*min(abs(x), abs(y));
    %z = log((1 + x*y)/(x + y));
end

%0.9375*