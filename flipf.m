function fnew = flipf(f) 
fnew = zeros(size(f));
for t = 1:length(f(1,1,1,:))
    for nE = 1:length(f(1,1,:,1))
        for nph = 1:length(f(:,1,1,1))
            fnew(nph,:,nE,t) = f(nph,end:-1:1,nE,t);
        end
        for nth = 1:length(f(1,:,1,1))
            fnew(:,nth,nE,t) = circshift(fnew(:,nth,nE,t),[16 0]);
        end
    end
end