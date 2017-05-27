function [xout, yout, zout]=rotsphere(x,y,z,Tht,Phi);

cT=cosd(Tht);
sT=sind(Tht);
cP=cosd(Phi);
sP=sind(Phi);


rot1=[[cT 0 sT];
   [0 1 0];
   [-sT 0 cT]];

rot2=[[cP -sP 0];
   [sP cP 0];
   [0 0 1]];


xo=rot1(1,1)*x+rot1(1,2)*y+rot1(1,3)*z;
yo=rot1(2,1)*x+rot1(2,2)*y+rot1(2,3)*z;
zo=rot1(3,1)*x+rot1(3,2)*y+rot1(3,3)*z;

xout=rot2(1,1)*xo+rot2(1,2)*yo+rot2(1,3)*zo;
yout=rot2(2,1)*xo+rot2(2,2)*yo+rot2(2,3)*zo;
zout=rot2(3,1)*xo+rot2(3,2)*yo+rot2(3,3)*zo;
