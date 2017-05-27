%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the perpendicular gradient in density from MMS
% observations. Inputs are the GSE Electric and Magnetic fields at the time
% slice (3 component vectors), the distribution function (evaluated at one 
% time, inputs are theta, phi, Energy), the vector of thetas, the vector of
% phis (both degrees), and the electron energy vector in eV. This assumes
% small drifts other than the ExB drift or near isotropy. 
% -BAW 10/23/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dn,varargout] = gradperpn(E,B,f,theta,phi,Enervec)
B = B*1e-9;
E = E*1e-3;
f = f*1e12;
e = 1.609e-19; m =9.11e-31;
b = log(Enervec(2))-log(Enervec(1));
Eperp = E - dot(E,B)/norm(B)^2*B;
[thetamat, phimat, Enermat] = ndgrid(theta,phi,e*Enervec);
Tmat = Enermat.*f;
phimat = phimat(:,:,1:end-1);
fpar = squeeze(mean(f(1,:,:),2))'+squeeze(mean(f(end,:,:),2))';
Weight = 2*pi^2*b*e*norm(B)/m^2/length(phi)/length(theta);
WE = b*sqrt(Enervec*e/m*2);
phiB = atan2(B(1),B(2));
xhat = [cos(phiB),-sin(phiB),0];
yhat = cross(B,xhat)/norm(B);
dn = -2*pi*e*Eperp/m*summ(fpar.*WE)...
    -2*Weight*(xhat*summ(logmean(Tmat(:,:,1:end-1),Tmat(:,:,2:end)).*sind(phimat))...
    -yhat*summ(logmean(Tmat(:,:,1:end-1),Tmat(:,:,2:end)).*cosd(phimat)));
F = sqrt(2/m^3)*squeeze(sum(sum(f.*sind(thetamat).*sqrt(Enermat),1),2))/length(phi)/length(theta)*2*pi^2;
n = diff(e*Enervec)*(F(1:end-1)+F(2:end))/2;
varargout{1} = n;
% dn = -2*pi*e*Eperp/m*summ(fpar.*WE)-2*(xhat*summ(f.*Wmat.*sind(phimat))-yhat*summ(f.*Wmat.*cosd(phimat)));
dn = dn*1e-8;
%dn in units of cm^-4;
