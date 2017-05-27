function [frot,varargout] = rotatef(f,thetain,phiin,thetaout,phiout,B)
%Here f should be squeezed down to f(phi,theta,E), all angles vectors
Ne = length(f(1,1,:));
frot = zeros(length(thetaout),length(phiout),Ne);
svar = size(thetaout);
if svar(1)>svar(2)
    phiout = phiout';
    thetaout = thetaout';
end
dt = thetain(2)-thetain(1);
dp = phiin(2)-phiin(1);
phivec = zeros(length(phiin)+2,1);
phivec(2:end-1) = phiin;
phivec(1) = phivec(2)-dp;phivec(end) = phivec(end-1)+dp;
fs = zeros(2+length(f(:,1,1)),length(f(1,:,1)));
[phioutm,thetaoutm] = meshgrid(phiout,thetaout);
[thetainterp,phiinterp] = Banglesinv(thetaoutm,phioutm,B);
[T,P] = meshgrid(double(thetain),double(phivec));
for i = 1:Ne
    fs(2:end-1,:) = squeeze(f(:,:,i));
    fs(end,:) = squeeze(f(1,:,i));
    fs(1,:) = squeeze(f(end,:,i));
    temp = interp2(T,P,fs,thetainterp,phiinterp);
    
    II = find(thetainterp<minn(thetain));
    if ~isempty(II)
        for ind = 1:length(II)
            thr = thetainterp(II(ind))/minn(thetain);
            weights = zeros(length(phiin),1);
            for j = 1:length(phiin)
                weights(j) = 1/(norm(thr*exp(1i*pi/180*phiinterp(ind))-exp(1i*pi/180*phiin(j)))+1e-6);
            end
            temp(II(ind)) = summ(weights.*squeeze(f(:,1,i)))./summ(weights);
        end
    end
    II = find(thetainterp>maxx(thetain));
    if ~isempty(II)
        for ind = 1:length(II)
            thr = (180-thetainterp(ind))/(180-maxx(thetain));
            weights = zeros(length(phiin),1);
            for j = 1:length(phiin)
                weights(j) = 1/(norm(thr*exp(1i*pi/180*phiinterp(II(ind)))-exp(1i*pi/180*phiin(j)))+1e-6);
            end
            temp(II(ind)) = summ(weights.*squeeze(f(:,end,i)))./summ(weights);
        end
    end
    
    frot(:,:,i)=temp;
end

varargout{1} = thetainterp;
varargout{2} = phiinterp;