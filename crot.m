function [phmat,thmat]=crot(ThtVec,PhiVec,Thtb,Phib)
thmat = zeros(length(ThtVec),length(PhiVec));
phmat = thmat;
for k=1:length(ThtVec)
    Tht=ThtVec(k);
    cT=cosd(Tht);
    sT=sind(Tht);

    z=cT+0*PhiVec;
    x=sT*cosd(PhiVec);
    y=sT*sind(PhiVec);

    [xout, yout, zout]=rotsphere(x,y,z,Thtb,Phib);
    [ph,th,~]=cart2sph(xout,yout,zout);
    th=th*180/pi;
    ph=ph*180/pi;

    phmat(k,:) = ph;
    thmat(k,:) = th;
end
thmat=90-thmat;

II=find(phmat<0);
phmat(II)=phmat(II)+360;


%xlabel('x')
%ylabel('y')
%zlabel('z')


%figure(4),clf
%for k=1:length(ThtVec)
%   plot(thmat(k,:),phmat(k,:),'g')
%   hold on
%end


%figure(4),hold on

%plot3(xout,yout,zout,'r')
%xlabel('x')
%ylabel('y')
%zlabel('z')