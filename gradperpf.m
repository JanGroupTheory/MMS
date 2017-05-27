%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the perpendicular gradient in density from MMS
% observations. Inputs are the GSE Electric and Magnetic fields at the time
% slice (3 component vectors), the distribution function (evaluated at one 
% time, inputs are theta, phi, Energy), the vector of thetas, the vector of
% phis (both degrees), and the electron energy vector in eV. This assumes
% small drifts other than the ExB drift or near isotropy. 
% -BAW 10/23/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [df,dfn,magE,maganf,ratioEf] = gradperpf(E,B,f,theta,phi,Enervec,Tperp)
B = B*1e-9;
E = E*1e-3;
f = f*1e12;
e = 1.609e-19; m =9.11e-31;
Tperp = Tperp*e;
Eperp = E - dot(E,B)/norm(B)^2*B;
[thetamat,phimat,Enermat] = ndgrid(theta,phi,e*Enervec);

rho = sqrt(2*m*Enermat)/e/norm(B).*(sind(thetamat)+1e-10);
phiB = atan2(B(1),B(2));
xhat = [cos(phiB),-sin(phiB),0];
yhat = cross(B,xhat)/norm(B);
maganf = 2*sqrt(squeeze(mean(f./rho.*sind(phimat),2).^2+mean(f./rho.*cosd(phimat),2).^2));
magE = e*norm(Eperp)/Tperp*squeeze(mean(f,2));
ratioEf = length(phi)/2*e*norm(Eperp)/Tperp*squeeze(mean(f,2))./sqrt(squeeze(sum(f./rho.*sind(phimat),2).^2+sum(f./rho.*cosd(phimat),2).^2));
dnx = -e*Eperp(1)/Tperp*squeeze(mean(f,2))-2/length(phi)*(xhat(1)*squeeze(sum(f./rho.*sind(phimat),2))-yhat(1)*squeeze(sum(f./rho.*cosd(phimat),2)));
dny = -e*Eperp(2)/Tperp*squeeze(mean(f,2))-2/length(phi)*(xhat(2)*squeeze(sum(f./rho.*sind(phimat),2))-yhat(2)*squeeze(sum(f./rho.*cosd(phimat),2)));
dnz = -e*Eperp(3)/Tperp*squeeze(mean(f,2))-2/length(phi)*(xhat(3)*squeeze(sum(f./rho.*sind(phimat),2))-yhat(3)*squeeze(sum(f./rho.*cosd(phimat),2)));
dn = zeros(3,50,length(Enervec));
dn(1,:,:) = dnx; dn(2,:,:) = dny; dn(3,:,:) = dnz;
df = dn*1e-8*1e-6;
normmat = zeros(size(df));
for i =1:3
    normmat(i,:,:)=squeeze(mean(f,2));
end
dfn = df./normmat*1e12;
%dn in units of cm^-4;
