function [theta0,phi0] = Banglesinv(theta,phi,B)
B = B/norm(B);
thetaB = acosd(B(3));
phiB = atan2d(B(1),B(2));
P1 = sind(theta).*cosd(phi);
P2 = sind(theta).*sind(phi);
Bh = cosd(theta);
Rmat = [cosd(phiB),sind(phiB),0;-sind(phiB),cosd(phiB),0;0,0,1]*...
    [1,0,0;0,cosd(thetaB),sind(thetaB);0,-sind(thetaB),cosd(thetaB)];
X = zeros(size(P1)); Y = zeros(size(P1)); Z = zeros(size(P1));
for i = 1:length(P1(:))
    temp = Rmat*[P1(i);P2(i);Bh(i)];
    X(i) = temp(1);
    Y(i) = temp(2);
    Z(i) = temp(3);
end
[phi0,theta0,~] = cart2sph(X,Y,Z);
phi0 = phi0*180/pi; theta0 = theta0*180/pi;
phi0(phi0<0) = phi0(phi0<0)+360;
theta0 = 90 - theta0;