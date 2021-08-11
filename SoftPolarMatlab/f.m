function z = f(x, y)
z = sign(x)*sign(y)*min(abs(x), abs(y));
end
