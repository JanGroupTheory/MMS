function getmmsdata(directory,tstart,tend,varargin) %tstart and tend should be seconds after the start of the file
DoPlot = 1;
if nargin>3
    if strcmp(varargin{1},'plot')
        DoPlot = 1;
    end
end
me=9.1e-31;
e=1.602e-19;
tstart = tstart/3600/24; tend = tend/3600/24;
tic
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];
for mms=1:4;
    Edataset = dir(['mms' num2str(mms) '*dce*']);
    fdataset=dir(['mms' num2str(mms) '*des-dist*']);
    mdataset=dir(['mms' num2str(mms) '*des-moms*']); 
    Bdataset=dir(['mms' num2str(mms) '*fgm*']);

    tn = spdfcdfread(mdataset.name,'Variables','Epoch');
    if mms == 1
        tinterp = tn;
    end
    try 
        temp = spdfcdfread(mdataset.name,'Variables',['mms' num2str(mms) '_des_numberdensity_dbcs_brst']);
    catch
        temp = spdfcdfread(mdataset.name,'Variables',['mms' num2str(mms) '_des_numberdensity_brst']);
    end
    n_esys(:,mms) = interp1(tn,temp,tinterp);
    RGSE = spdfcdfread(Bdataset.name,'Variables',['mms' num2str(mms) '_fgm_r_gse_brst_l2']);
    tBin = spdfcdfread(Bdataset.name,'Variables','Epoch');
    if mms == 1
        tbint = tBin;
    end
    temp = spdfcdfread(Bdataset.name,'Variables',['mms' num2str(mms) '_fgm_b_gse_brst_l2']);
    for di = 1:4
        BGSEsys(:,di,mms) = interp1(tBin,temp(:,di),tinterp);
    end
    tR = spdfcdfread(Bdataset.name,'Variables','Epoch_state');
    X(:,mms) = interp1(tR,RGSE(:,1),tinterp);
    Y(:,mms) = interp1(tR,RGSE(:,2),tinterp);
    Z(:,mms) = interp1(tR,RGSE(:,3),tinterp);
    
    tvecf=spdfcdfread(fdataset.name,'Variables','Epoch');
    if mms == 1
        t0=tvecf(1);
        maxt = max(tvecf-t0);
    end
    tvecf=tvecf-t0;
    if tend ==0
        tend = maxt;
    end
%     fdata=spdfcdfread(fdataset.name);
%     Bdata=spdfcdfread(Bdataset.name);
    tvecb=spdfcdfread(Bdataset.name,'Variables','Epoch');
    tvecb=tvecb-t0;
    B=spdfcdfread(Bdataset.name,'Variables',['mms' num2str(mms) '_fgm_b_gse_brst_l2']);
    xyz=spdfcdfread(Bdataset.name,'Variables',['mms' num2str(mms) '_fgm_r_gse_brst_l2']);
    xyz=xyz(1,:);
    f=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_dist_brst']);
    f = flipf(f);
    
    fphi=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_phi_brst']);  % phi angles of distribution data (in GSE)
    ftht=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_theta_brst']); % theta angles of distribution data (in GSE)
    try
        Ener2=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_energy0_brst']);
        Ener1=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_energy1_brst']);
    catch
        Ener2=spdfcdfread(fdataset.name,'Variables',['mms' num2str(mms) '_des_energy_brst']);
        Ener1=Ener2;
    end

    E = spdfcdfread(Edataset.name,'Variables',['mms' num2str(mms) '_edp_dce_gse_brst_l2']);
    tE = spdfcdfread(Edataset.name,'Variables',['mms' num2str(mms) '_edp_epoch_brst_l2']);
    tE = tE - t0;

    E = E(1:10:end,:); %Avoid same time in interpolation
    tE = tE(1:10:end);

    Bf=zeros([length(tvecf),4]); % prepare to interpolate B to same time vector as distribution functions. 
    Ef = Bf;
    if tvecb(1)>tvecf(1)
        tvecb(1)=tvecf(1);
    end
    if tvecb(end)<tvecf(end)
        tvecb(end)=tvecf(end);
    end
    if tE(1)>tvecf(1)
        tE(1)=tvecf(1);
    end
    if tE(end)<tvecf(end)
        tE(end)=tvecf(end);
    end
    for k=1:4
        Bf(:,k)=interp1(tvecb,B(:,k),tvecf);% interpolate B to same time vector as distribution functions.
        if k < 4
            Ef(:,k)=interp1(tE,E(:,k),tvecf);% interpolate E to same time vector as distribution functions.
        else
            Ef(:,4) = sqrt(Ef(:,1).^2 + Ef(:,2).^2 + Ef(:,3).^2);
        end
    end
    
    Tepar = spdfcdfread(mdataset.name,'Variables',['mms' num2str(mms) '_des_temppara_brst']);
    Teperp = spdfcdfread(mdataset.name,'Variables',['mms' num2str(mms) '_des_tempperp_brst']);
    Te = (Tepar+2*Teperp)/3;
    
    ThtVec=linspace(0,1,50)*180;   %field aligned coords
    PhiVec=linspace(0,2,101)*180;
    PhiVec=PhiVec(1:100);

    Ntvec = find(((tvecf<tend).*(tvecf>tstart))==1)';

    IIe=1:17; % index of energies to be used

    EnerVec=Ener2(IIe);
    vz1 = zeros(length(IIe),length(ThtVec));
    vp1 = vz1; vz2 = vz1; vp2 = vp1;
    for kk=IIe
        vz1(kk,:)=sqrt(EnerVec(kk))*cos(ThtVec*pi/180); %parallel velocity for even timesteps
        vp1(kk,:)=sqrt(EnerVec(kk))*sin(ThtVec*pi/180);  %perp velocity for even timesteps
    end
    EnerVec=Ener1(IIe);
    for kk=IIe
        vz2(kk,:)=sqrt(EnerVec(kk))*cos(ThtVec*pi/180); %parallel velocity for odd timesteps
        vp2(kk,:)=sqrt(EnerVec(kk))*sin(ThtVec*pi/180); %perp velocity for odd timesteps
    end

    vz1=sqrt(2*e/me)*vz1/1e7;
    vp1=sqrt(2*e/me)*vp1/1e7;
    vz2=sqrt(2*e/me)*vz2/1e7;
    vp2=sqrt(2*e/me)*vp2/1e7;

    Cphi=zeros([length(ThtVec), length(PhiVec),length(IIe)]); % matrices useful for calculating agyrotropy
    Sphi=zeros([length(ThtVec), length(PhiVec),length(IIe)]);
    for kk=1:length(ThtVec)
        for ll=1:length(IIe)
            Sphi(kk,:,ll)=sin(PhiVec*pi/180);
            Cphi(kk,:,ll)=cos(PhiVec*pi/180);
        end
    end


    Fave=zeros([length(ThtVec), length(IIe),length(Ntvec)]);  % log10 of gyro averaged f
    FaveS=zeros([length(ThtVec), length(IIe),length(Ntvec)]); % sin-dependency of agyrotropy
    FaveC=zeros([length(ThtVec), length(IIe),length(Ntvec)]); % cosine-dependency of agyrotropy
    phmid=zeros([100 length(Ntvec)]); % phi in GSE for points perp to B
    thmid=zeros([100 length(Ntvec)]); % theta in GSE for points perp to B
    phB=zeros(length(Ntvec)); % phi in GSE of  B
    thB=zeros(length(Ntvec)); % theta in GSE of B
    gradf = zeros(length(tvecf),3,50,length(IIe));
    gradfn = gradf;
    ratioEf = zeros(length(tvecf),50,length(IIe));
    magEf = ratioEf; maganf = magEf;
    fphiav = zeros(length(Ntvec),length(ThtVec),length(IIe));
    gradn = zeros(length(tvecf),3);
%     ue = gradn;
    cnt=0;
    for Nft=Ntvec
        cnt=cnt+1;

        Bxyz=double(Bf(Nft,1:3));
%         Bu=Bxyz/sqrt(Bxyz*Bxyz');
%         PhiB=atan2d(Bu(2),Bu(1)); % Phi angle of magnetic field in GSE coords
%         if PhiB < 0
%             PhiB=PhiB+360;
%         end
%         ThtB=acosd(Bu(3));  % Theta angle of magnetic field in GSE coords
        farr=double(squeeze(f(:,:,IIe,Nft)));  % select data at time Nft
        fphiplt=squeeze(fphi(Nft,:));


%         [phmat,thmat]=crot(ThtVec,PhiVec,ThtB,PhiB);  %calculate points for interpolating f into field alinged coords. 
%         thmid(:,cnt)=squeeze(thmat(25,:))';
%         phmid(:,cnt)=squeeze(phmat(25,:))';
%         thB(cnt)=ThtB;
%         phB(cnt)=PhiB;


%         Ner=length(IIe);
%         fmat=rotatef(thmat,phmat,ftht,fphiplt,farr,Ner);
        fmat=rotatef(farr,ftht,fphiplt,ThtVec,PhiVec,Bxyz);

        FmatLog=log10(fmat+1e-31);
        Fplt=squeeze(mean(FmatLog(:,:,:),2));
        Fsplt=squeeze(mean(Sphi.*FmatLog(:,:,:),2));
        Fcplt=squeeze(mean(Cphi.*FmatLog(:,:,:),2));
        
%         if DoPlot
%             figure(4),clf
%             IIp=1:14;
%             subplot(2,1,1);
%             if mod(Nft,2)==0
%             pcolor(vz1',vp1',Fplt+sqrt(Fsplt.^2+Fcplt.^2)),
%             hold on
%             pcolor(vz1',-vp1',Fplt-sqrt(Fsplt.^2+Fcplt.^2)),
%             else
%                 pcolor(vz2',vp2',Fplt+sqrt(Fsplt.^2+Fcplt.^2)),
%                 hold on
%                 pcolor(vz2',-vp2',Fplt-sqrt(Fsplt.^2+Fcplt.^2)),
%             end
%             shading interp,
%             caxis([-28.6 -25])
%             axis([-1 1 0 1]* vz1(IIp(end),1))
% 
% 
%             subplot(2,1,2)
%             if mod(Nft,2)==0
%                 pcolor(vz1(IIp,:)',vp1(IIp,:)',smoo(smoo(sqrt(Fsplt(:,IIp).^2+Fcplt(:,IIp).^2)',2)',2)),
%             else
%                 pcolor(vz2(IIp,:)',vp2(IIp,:)',smoo(smoo(sqrt(Fsplt(:,IIp).^2+Fcplt(:,IIp).^2)',2)',2)),
%             end
%             shading interp
%             caxis([0 0.35])
%             axis([-1 1 0 1]* vz1(IIp(end),1))
%             saveas(gcf,[directory '\Dists\MMS' num2str(mms) '_' num2str(Nft) 't.png']);
%         end

        Fave(:,:,cnt)=squeeze(mean(FmatLog(:,:,:),2));
        FaveS(:,:,cnt)=squeeze(mean(Sphi.*FmatLog(:,:,:),2));
        FaveC(:,:,cnt)=squeeze(mean(Cphi.*FmatLog(:,:,:),2));

        drawnow
        
        if mod(Nft,2) == 0
            Enervec = Ener2;
        else
            Enervec = Ener1;
        end
        
        fphiav(Nft,:,:) = squeeze(mean(fmat,2));
        [gradf(Nft,:,:,:),gradfn(Nft,:,:,:),magEf(Nft,:,:),maganf(Nft,:,:),ratioEf(Nft,:,:)] = gradperpf(Ef(Nft,1:3),Bf(Nft,1:3),fmat,ThtVec,PhiVec,Enervec(IIe),Te(Nft));
        gradn(Nft,:) = gradperpn(Ef(Nft,1:3),Bf(Nft,1:3),fmat(:,:,5:17),ThtVec,PhiVec,Enervec(5:17));
%         ue(Nft,:) = intue(farr,ftht,fphiplt,Enervec(IIe));
    end
    
    %f gradient
    perpplanegrad = squeeze(mean(gradf(:,:,23:27,:),3));
    perpplanegradn = squeeze(mean(gradfn(:,:,23:27,:),3));
    perpplaneratio = squeeze(mean(ratioEf(:,23:27,:),2));
    ppE = squeeze(mean(magEf(:,23:27,:),2));
    ppan = squeeze(mean(maganf(:,23:27,:),2));
    AvEPerpGrad = squeeze(mean(perpplanegrad(:,:,9:13),3));
    AvErat = squeeze(mean(perpplaneratio(:,9:13),2));
    AvppE = squeeze(mean(ppE(:,9:13),2));
    Avppan = squeeze(mean(ppan(:,9:13),2));
    gradmag = sqrt(AvEPerpGrad(:,1).^2+AvEPerpGrad(:,2).^2+AvEPerpGrad(:,3).^2);
    
    toc
    tic
    varlist = {'Bf','Ef','Ener1','Ener2','EnerVec','Fave','FaveC','FaveS',...
        'fphi','ftht','Ntvec','PhiVec','t0','tvecb','tE','tvecf',...
        'ThtVec','vp1','vp2','vz1','vz2','gradn','gradf','gradmag',...
        'perpplanegrad','perpplaneratio','ppE','ppan','AvEPerpGrad','AvErat',...
        'AvppE','Avppan','gradfn','perpplanegradn'};
    save([Datadir '\MMS' num2str(mms) 'ave'],varlist{:});
end

for nt = 1:length(n_esys(:,1))
    B = squeeze(mean(BGSEsys(nt,1:3,:),3));
    magB = sqrt(summ(B.*B));
    B = B/magB;
    mmsdens = n_esys(nt,:);
    mmsx = X(nt,:)-meann(X(nt,:));
    mmsy = Y(nt,:)-meann(Y(nt,:));
    mmsz = Z(nt,:)-meann(Z(nt,:));
    xs = summ(mmsx.^2);
    ys = summ(mmsy.^2);
    zs = summ(mmsz.^2);
    xy = summ(mmsx.*mmsy);
    xz = summ(mmsx.*mmsz);
    yz = summ(mmsz.*mmsy);
    A = [xs xy xz; xy ys yz; xz yz zs];
    b = [summ(mmsx.*mmsdens);summ(mmsy.*mmsdens);summ(mmsz.*mmsdens)];
    Gradn(nt,:) = A\b;
%     altdn(nt,1) = sum((mmsdens-meann(mmsdens)).*mmsx)./(zs);
%     altdn(nt,2) = sum((mmsdens-meann(mmsdens)).*mmsy)./(ys);
%     altdn(nt,3) = sum((mmsdens-meann(mmsdens)).*mmsz)./(zs);
    Gradperpn(nt,:) = Gradn(nt,:) - dot(Gradn(nt,:) ,B)*B; 
end
 save([Datadir '\MMS_FD'],'Gradn','Gradperpn','n_esys','BGSEsys');