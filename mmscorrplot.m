function mmscorrplot(input,varargin) %input should be matrix of 4 column vectors of equal length
% input = input - meann(input);
if nargin>1
    npoints = varargin{1};
else
    npoints = 10;
end
[cr,lgs]=xcorr(input,npoints,'unbiased');
crp = zeros(size(cr));
for row = 1:4
    for col = 1:4
        %renormalization parameters
        Ar = sqrt(maxx(cr(:,4*(row-1)+row)));
        Ac = sqrt(maxx(cr(:,4*(col-1)+col)));
        nm = 4*(row-1)+col;
        crp(:,nm)=cr(:,nm)/(Ar*Ac);
        subplot(4,4,nm)
        stem(lgs,crp(:,nm),'.')
        title(sprintf('c_{%d%d}',row,col))
        ylim([minn(crp(:,nm)) 1])
    end
end