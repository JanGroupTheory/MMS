function mmsplotlayer(directory, v, Xwidth, varargin)
havetimes=false;
DoPlot = false;
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];
Bave=0;
xyzMat=zeros([3,4]);
%width of event is about 100 km
Nsb=900; %start sample consider
Nse=1500; %start sample consider
Nskip=5;
IIe = 9;

datafile = dir([Datadir '\*' num2str(1) '*lmn*']);
datafile = [Datadir '\' datafile.name];
load(datafile);
datafile = dir([Datadir '\*' num2str(1) '*ave*']);
datafile = [Datadir '\' datafile.name];
load(datafile);

mastert = tvecf;
Bfln = Bfl;

Bxbar = zeros(length(Bf(:,1)),1);Bx = zeros(length(Bf(:,1)),4);
By=Bx;Bz=By;vx=zeros(length(Blmn(:,1)),4);vy=vx;vz=vy;
nemat = Bx; vxf=Bx;vyf=vxf;vzf=vyf;Pavef=nemat; 
for mms = 1:4
    if mms>1
        datafile = dir([Datadir '\*' num2str(mms) '*lmn*']);
        datafile = [Datadir '\' datafile.name];
        load(datafile);
        datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
        datafile = [Datadir '\' datafile.name];
        load(datafile);
        for j = 1:3
            Bfln(:,j) = interp1(tvecf,Bfl(:,j),mastert);
        end
        Bfl = Bfln;
    end
    Bave = Bave+Bfl/4;
    Bxbar = Bxbar + Bfl(:,1)/4;
    Bx(:,mms) = Bfl(:,1);
    By(:,mms) = Bfl(:,2);
    Bz(:,mms) = Bfl(:,3);
    nemat(:,mms) = interp1(tB,ne(:),mastert); %(:) avoids confusion with function
    vxf(:,mms) = interp1(tB,vlmn(:,1),mastert);
    vyf(:,mms) = interp1(tB,vlmn(:,2),mastert);
    vzf(:,mms) = interp1(tB,vlmn(:,3),mastert);
    Pavef(:,mms) =1/3*interp1(tB,Ppar+Pperp1+Pperp2,mastert); 
    xyzMat(:,mms)=mean(rlmn(:,1:3),1);
end

II = Nsb:Nse;
Bave = Bave(II,:);
Bx = Bx(II,:);
By = By(II,:);
Bz = Bz(II,:);
vx = vxf(II,:);
vy = vyf(II,:);
vz = vzf(II,:);
tf = mastert(II);
nemat = nemat(II,:);
vz = vz + v;
figure()
plot(secssince(tf),Bx,'b',secssince(tf),By,'r',secssince(tf),Bz,'g',secssince(tf),sqrt(Bx.^2+By.^2+Bz.^2),'k')
legend('Bx','By','Bz','|B|');
title('B vs t');


for k=1:3
    xyzMat(k,:)=xyzMat(k,:)-mean(xyzMat(k,:)); %relative xyz of the four spacecraft
end


Bxmax=maxx(abs(Bx));
ds=30e-3*v; %distance per sample
% ds=1;
Zcm = -Bave(:,1)'/Bxmax*Xwidth;
% Xcm=smoo(Xcm,30)';
Xcm=(Nsb:Nse)*ds;
XXcm=[Xcm' Xcm'];
ZZcm=[Zcm' Zcm'+2];

if nargin>3
    for index = 1:length(varargin)
        if strcmp(varargin{index},'classic')
            Nts = round(varargin{index+1})-Nsb+1;
            havetimes = true;
        end
    end
end

figure(5),clf
subplot(2,1,1)
col='brgk';
for k=1:4
    plot(Xcm+xyzMat(1,k),Zcm+xyzMat(3,k),col(k))
    hold on
    axis equal
    poss = [Nsb Nse]*ds;
    xlim([minn(poss) maxx(poss)]);
end
subplot(2,1,2)
for k=1:4
    pcolor(XXcm+xyzMat(1,k),ZZcm+xyzMat(3,k),[By(:,k) By(:,k)])
    shading interp
    hold on
    quiver(Xcm+xyzMat(1,k),Zcm+xyzMat(3,k),Bx(:,k)',Bz(:,k)',col(k))
    if havetimes
        plot((Xcm(Nts(k,:))+xyzMat(1,k)),Zcm(Nts(k,:))+xyzMat(3,k) ,['o' col(k)],'markersize',10,'linewidth',2.5)
    end  
    axis equal
    xlim([minn(poss) maxx(poss)]);
    colorbar
end
if nargin>3
    if DoPlot
        figure()
        for mms = 1:4
            gradnp = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);

            for k = 1:3
                gradnp(:,k) = interp1(tvecf,gradnl(:,k),mastert);
            end
            gradnp = gradnp(II,:);
            for q = 1:length(gradnp(:,1))
                gradnpn(q,:) = gradnp(q,:)/norm(gradnp(q,:));
            end
            subplot(2,1,1)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[gradnp(:,2) gradnp(:,2)])
            shading interp
            caxis([-1 1]*maxx(sqrt(gradnp(:,1).^2+gradnp(:,2).^2+gradnp(:,3).^2)))
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),gradnp(1:Nskip:end,1)',gradnp(1:Nskip:end,3)',col(mms))
            axis equal
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            xlim([minn(poss) maxx(poss)]);
            colorbar
            subplot(2,1,2)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[gradnpn(:,2) gradnpn(:,2)])
            shading interp
            caxis([-1 1])
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),gradnpn(1:Nskip:end,1)',gradnpn(1:Nskip:end,3)',col(mms))
            axis equal
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            xlim([minn(poss) maxx(poss)]);
            colorbar
        end
    end
    DoPlot = false;
    for index = 1:length(varargin)
        if strcmp(varargin{index},'gradf')||strcmp(varargin{index},'all')
            DoPlot = true;
        end
    end
    if DoPlot
        figure()
        for mms = 1:4
            gradnp = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);

            for k = 1:3
                gradnp(:,k) = interp1(tvecf,squeeze(mean(perpplanegradl(:,k,IIe),3)),mastert);
            end
            gradnp = gradnp(II,:);
            for q = 1:length(gradnp(:,1))
                gradnpn(q,:) = gradnp(q,:)/norm(gradnp(q,:));
            end
            subplot(2,1,1)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[gradnp(:,2) gradnp(:,2)])
            shading interp
            caxis([-1 1]*maxx(sqrt(gradnp(:,1).^2+gradnp(:,2).^2+gradnp(:,3).^2)))
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),gradnp(1:Nskip:end,1)',gradnp(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
            subplot(2,1,2)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[gradnpn(:,2) gradnpn(:,2)])
            caxis([-1 1])
            shading interp
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),gradnpn(1:Nskip:end,1)',gradnpn(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
        end
    end
    DoPlot = false;
    for index = 1:length(varargin)
        if strcmp(varargin{index},'anis')||strcmp(varargin{index},'all')
            DoPlot = true;
        end
    end
    if DoPlot
        for mms = 1:4
            figure(43)
            Nskip = 20;
            gradnp = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);

            Fz = smoo(squeeze(mean(Fave(1:5,IIe,:),1)),2);                
            Fp = smoo(squeeze(mean(Fave(23:27,IIe,:),1)),2);  
            Anis = interp1(tvecf(Ntvec),(Fz-Fp),mastert);
            Anis = Anis(II);
            for k = 1:3
                gradnp(:,k) = interp1(tvecf,squeeze(mean(perpplanegradl(:,k,IIe),3)),mastert);
            end
            gradnp = gradnp(II,:);
            if mms == 1
                Anisscale = (maxx(XXcm)-minn(XXcm))/(maxx(Anis)-minn(Anis))/4;
            end
%                 tubeplot(Xcm'+xyzMat(1,mms),Zcm'+xyzMat(3,mms),Anisscale*Anis,1,Anis)
            tubeplot(Xcm'+xyzMat(1,mms),Zcm'+xyzMat(3,mms),Anisscale*Anis,1,mms*ones(size(Anis)))
%                 surf(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),Anisscale*[Anis Anis+.02*maxx(abs(Anis))])
            shading interp
            caxis([0 5])
%                 caxis([-1 1]*maxx(abs(Anis)))
            hold on
%                 quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),gradnp(1:Nskip:end,1)',gradnp(1:Nskip:end,3)',col(mms))
            if havetimes
                plot3((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,Anisscale*Anis(Nts(mms,:)),['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
            figure(44)
            subplot(2,1,1)
            plot(Xcm'+xyzMat(1,mms),Zcm'+xyzMat(3,mms))
            hold on
            subplot(2,1,2)
            plot(Xcm'+xyzMat(1,mms),Anisscale*Anis)
            hold on
        end
    end
    DoPlot = false;
    for index = 1:length(varargin)
        if strcmp(varargin{index},'vi')||strcmp(varargin{index},'all')
            DoPlot = true;
        end
    end
    if DoPlot
        for mms = 1:4
            figure(45)
            Nskip = 20;
            vlmns = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);
            for i = 1:3
                vlmns(:,i) = interp1(tB-t0,vlmn(:,i),mastert);
            end
            vlmns = vlmns(II,:);
            vlmns(:,1)=vlmns(:,1)+v;
            vlmns(:,2) = vlmns(:,2) - meann(vlmns(:,2));
            vlmnsn = vlmns.*[(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2)];
            subplot(2,1,1)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmns(:,2) vlmns(:,2)])
            title('u_i');
            shading interp
            caxis([-1 1]*maxx(sqrt(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2)))
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmns(1:Nskip:end,1)',vlmns(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
            subplot(2,1,2)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmnsn(:,2) vlmnsn(:,2)])
            shading interp
            caxis([-1 1])
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmnsn(1:Nskip:end,1)',vlmnsn(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
        end
    end
    DoPlot = false;
    for index = 1:length(varargin)
        if strcmp(varargin{index},'J')||strcmp(varargin{index},'all')
            DoPlot = true;
        end
    end
    if DoPlot
        for mms = 1:4
            figure(46)
            Nskip = 5;
            vlmns = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);
            for i = 1:3
                vlmns(:,i) = interp1(tB-t0,Jlmn(:,i),mastert);
            end
            vlmns = vlmns(II,:);
            vlmnsn = vlmns.*[(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2)];
            subplot(2,1,1)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmns(:,2) vlmns(:,2)])
            title('J');
            shading interp
            caxis([-1 1]*maxx(sqrt(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2)))
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmns(1:Nskip:end,1)',vlmns(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
            subplot(2,1,2)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmnsn(:,2) vlmnsn(:,2)])
            shading interp
            caxis([-1 1])
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmnsn(1:Nskip:end,1)',vlmnsn(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
        end
    end
    DoPlot = false;
    for index = 1:length(varargin)
        if strcmp(varargin{index},'ue')||strcmp(varargin{index},'all')
            DoPlot = true;
        end
    end
    if DoPlot
        for mms = 1:4
            figure(47)
            Nskip = 5;
            vlmns = zeros(length(mastert),3);
            datafile = dir([Datadir '\*' num2str(mms) '*ave*']);
            datafile = [Datadir '\' datafile.name];
            load(datafile);
            for i = 1:3
                vlmns(:,i) = interp1(tB-t0,uelmn(:,i),mastert);
            end
            vlmns = vlmns(II,:);
            vlmnsn = vlmns.*[(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2),...
                (vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2).^(-1/2)];
            subplot(2,1,1)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmns(:,2) vlmns(:,2)])
            title('u_e');
            shading interp
            caxis([-1 1]*maxx(sqrt(vlmns(:,1).^2+vlmns(:,2).^2+vlmns(:,3).^2)))
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmns(1:Nskip:end,1)',vlmns(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
            subplot(2,1,2)
            pcolor(XXcm+xyzMat(1,mms),ZZcm+xyzMat(3,mms),[vlmnsn(:,2) vlmnsn(:,2)])
            shading interp
            caxis([-1 1])
            hold on
            quiver(Xcm(1:Nskip:end)+xyzMat(1,mms),Zcm(1:Nskip:end)+xyzMat(3,mms),vlmnsn(1:Nskip:end,1)',vlmnsn(1:Nskip:end,3)',col(mms))
            if havetimes
                plot((Xcm(Nts(mms,:))+xyzMat(1,mms)),Zcm(Nts(mms,:))+xyzMat(3,mms) ,['o' col(mms)],'markersize',10,'linewidth',2.5)
            end  
            axis equal
            xlim([minn(poss) maxx(poss)]);
            colorbar
        end
    end 
end
