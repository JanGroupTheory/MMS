%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Getting the [l m n] using min variance analysis
%%%%%%%%%% V=[time,Vx,Vy,Vz,Vl,Vm,Vn]
%%%%%%%%%% Blmn=[time,posx,posy,posz,Bl,Bm,Bn]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lmn,Blmn]=GSM_mva(B,v,yr,mn,dy,hrstart,mstart,hrend,mend)

t = [yr,mn,dy,hrstart,mstart,00];
dtnumi=datenum(t);
t = [yr,mn,dy,hrend,mend,00];
dtnumf=datenum(t);
        
temp = abs(B(:,1) - dtnumi);
[~,tBi]= min(temp);

temp = abs(B(:,1) - dtnumf);
[~,tBf]= min(temp);

temp = abs(v(:,1) - dtnumi);
[~,tvi]= min(temp);

temp = abs(v(:,1) - dtnumf);
[~,tvf]= min(temp);

%% Minimum Variance Analysis
va=v(tvi:tvf,3:5);
[a,~]=size(va);
I=ones(1,a);
va=I*va;
l=transpose(va./norm(va));

MB=zeros(3,3);
B1=B(tBi:tBf,5:7);
[a,b]=size(B1);
mean=zeros(1,b);
for i=1:a
mean=B1(i,:)+mean;
end
mean=mean/a;
meanp=[];
for k=1:a
    meanp=[meanp;mean]; 
end
MB=1/(a-1)*(B1-meanp)'*(B1-meanp);

[VB,DB]=eig(MB);
% disp(DB);
DB(~DB)=inf;
[~,min_loc] = min(min(DB));
[~,max_loc] = max(min(DB));
% disp(min_loc);
% disp(max_loc);
% med_loc = 6-min_loc-max_loc;
np = VB(:,max_loc);
n = cross(np,l);
if n(1)<0
    n = -n; 
end

z = [0;0;1];
l = (z - dot(z,n)*n)/norm(z - dot(z,n)*n);
m=cross(n,l);
lmn=[l m n];
%% Getting Blmn
 [s,~] = size(B);
 for i=1:s
 pos=B(i,2:4);
 Bxyz=B(i,5:7);
 Blmn_1=Bxyz*lmn;
 pos_lmn=pos*lmn;
 Blmn(i,2:7)=[pos_lmn(1,:),Blmn_1(1,:)];
 end
 Blmn(:,8)=B(:,8);
 Blmn(:,1)=B(:,1);

end
