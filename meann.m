function xbar = meann(x)
xbar = x;
while length(xbar)>1
    xbar = mean(xbar);
end