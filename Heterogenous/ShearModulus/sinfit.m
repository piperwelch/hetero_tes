% define function to fit to
function y = sinfit(x, a, dc, th0)
    y = a*sin(4*x + th0) + dc;
end