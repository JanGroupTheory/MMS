function a = summ(x)
a = x;
while max(size(a))>1
    a = squeeze(sum(a));
end