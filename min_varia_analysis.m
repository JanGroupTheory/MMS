%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Getting the [l m n] using min variance analysis
%%%%%%%%%% V=[time,Vx,Vy,Vz,Vl,Vm,Vn]
%%%%%%%%%% Blmn=[time,posx,posy,posz,Bl,Bm,Bn]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lmn,Blmn,V]=min_varia_analysis(B,v,yr,mn,dy,hrstart,mstart,hrend,mend)

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
disp(DB);
DB(~DB)=inf;
[~,min_loc] = min(min(DB));
[~,max_loc] = max(max(DB));
med_loc = 6-min_loc-max_loc;
n = VB(:,min_loc);
if n(3)<0
    n = -n; 
end

va=v(tvi:tvf,3:5);
[a,~]=size(va);
I=ones(1,a);
va=I*va;
l=transpose(va./norm(va));
l = (l - dot(l,n)*n)/norm(l - dot(l,n)*n);
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
%% Getting Vlmn
[s,~] = size(v);
 for i=1:s
 vxyz=v(i,3:5);
 vlmn_1=vxyz*lmn;
 V(i,2:7)=[vxyz,vlmn_1];
 end
 V(:,1)=v(:,1);
end
