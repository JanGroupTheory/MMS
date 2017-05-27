function mmsloadGSE(directory)
tic
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];
for i = 1:4
    savefile = [Datadir '\MMS' num2str(i) '_GSE'];
    vfile = dir(['mms' num2str(i) '*dis-moms*']);
    vfile = vfile.name;
    Pfile = dir(['mms' num2str(i) '*des-moms*']);
    Pfile = Pfile.name;
    filename = dir(['mms' num2str(i) '*fgm*']);
    filename = filename.name;
    ffile = dir(['mms' num2str(i) '*des-dist*']);
    ffile = ffile.name;
    tf =spdfcdfread(ffile,'Variables','Epoch');
    Efile = dir(['mms' num2str(i) '*edp*dce*']);
    Efile = Efile.name;
    Evar = ['mms' num2str(i) '_edp_dce_gse_brst_l2'];
    tEvar = ['mms' num2str(i) '_edp_epoch_brst_l2'];
    rvar = ['mms' num2str(i) '_fgm_r_gse_brst_l2'];
    nvar = ['mms' num2str(i) '_des_numberdensity_dbcs_brst'];
    nivar = ['mms' num2str(i) '_dis_numberdensity_dbcs_brst'];
    Bvar = ['mms' num2str(i) '_fgm_b_gse_brst_l2'];
    vvar = ['mms' num2str(i) '_dis_bulkspeed_dbcs_brst'];
    vvarx = ['mms' num2str(i) '_dis_bulkx_dbcs_brst'];
    vvary = ['mms' num2str(i) '_dis_bulky_dbcs_brst'];
    vvarz = ['mms' num2str(i) '_dis_bulkz_dbcs_brst'];
    uevar = ['mms' num2str(i) '_des_bulkspeed_dbcs_brst'];
    uevarx = ['mms' num2str(i) '_des_bulkx_dbcs_brst'];
    uevary = ['mms' num2str(i) '_des_bulky_dbcs_brst'];
    uevarz = ['mms' num2str(i) '_des_bulkz_dbcs_brst'];
    Pvarxx = ['mms' num2str(i) '_des_presxx_dbcs_brst'];
    Pvarxy = ['mms' num2str(i) '_des_presxy_dbcs_brst'];
    Pvarxz = ['mms' num2str(i) '_des_presxz_dbcs_brst'];
    Pvaryy = ['mms' num2str(i) '_des_presyy_dbcs_brst'];
    Pvaryz = ['mms' num2str(i) '_des_presyz_dbcs_brst'];
    Pvarzz = ['mms' num2str(i) '_des_preszz_dbcs_brst'];
    Tparvar = ['mms' num2str(i) '_des_temppara_brst'];
    Tperpvar = ['mms' num2str(i) '_des_tempperp_brst'];
    Tiparvar = ['mms' num2str(i) '_dis_temppara_brst'];
    Tiperpvar = ['mms' num2str(i) '_dis_tempperp_brst'];
    
    rGSEslow = double(spdfcdfread(filename,'Variables',rvar));
    tr = spdfcdfread(filename,'Variables','Epoch_state');
    BGSE = double(spdfcdfread(filename,'Variables',Bvar));
    EGSEf = double(spdfcdfread(Efile,'Variables',Evar));
    try
        vGSEslow = double(spdfcdfread(vfile,'Variables',vvar));
        vGSEslowx = double(spdfcdfread(vfile,'Variables',vvarx));
        vGSEslowy = double(spdfcdfread(vfile,'Variables',vvary));
        vGSEslowz = double(spdfcdfread(vfile,'Variables',vvarz));
        vGSEslow = [vGSEslowx,vGSEslowy,vGSEslowz,vGSEslow];
    catch
        vvar = ['mms' num2str(i) '_dis_bulkv_dbcs_brst'];
        vGSEslow = double(spdfcdfread(vfile,'Variables',vvar)); 
    end

    try
        ueGSEslow = double(spdfcdfread(Pfile,'Variables',uevar));
        ueGSEslowx = double(spdfcdfread(Pfile,'Variables',uevarx));
        ueGSEslowy = double(spdfcdfread(Pfile,'Variables',uevary));
        ueGSEslowz = double(spdfcdfread(Pfile,'Variables',uevarz));
        ueGSEslow = [ueGSEslowx,ueGSEslowy,ueGSEslowz,ueGSEslow];
    catch
        uevar = ['mms' num2str(i) '_des_bulkv_dbcs_brst'];
        ueGSEslow = double(spdfcdfread(Pfile,'Variables',uevar));
    end
    try
        Pxxs = double(spdfcdfread(Pfile,'Variables',Pvarxx));
        Pxys = double(spdfcdfread(Pfile,'Variables',Pvarxy));
        Pxzs = double(spdfcdfread(Pfile,'Variables',Pvarxz));
        Pyys = double(spdfcdfread(Pfile,'Variables',Pvaryy));
        Pyzs = double(spdfcdfread(Pfile,'Variables',Pvaryz));
        Pzzs = double(spdfcdfread(Pfile,'Variables',Pvarzz));
    catch
        Ptens = double(spdfcdfread(Pfile,'Variables',['mms' num2str(i) '_des_prestensor_dbcs_brst']));
        Pxxs = squeeze(Ptens(1,1,:));
        Pxys = squeeze(Ptens(1,2,:));
        Pxzs = squeeze(Ptens(1,3,:));
        Pyys = squeeze(Ptens(2,2,:));
        Pyzs = squeeze(Ptens(2,3,:));
        Pzzs = squeeze(Ptens(3,3,:));  
    end
    Tpars = double(spdfcdfread(Pfile,'Variables',Tparvar));
    Tperps = double(spdfcdfread(Pfile,'Variables',Tperpvar));
    Tipars = double(spdfcdfread(vfile,'Variables',Tiparvar));
    Tiperps = double(spdfcdfread(vfile,'Variables',Tiperpvar));
    try
        nes = double(spdfcdfread(Pfile,'Variables',nvar));
    catch
        nvar = ['mms' num2str(i) '_des_numberdensity_brst'];
        nes = double(spdfcdfread(Pfile,'Variables',nvar));
    end
    try
        nis = double(spdfcdfread(vfile,'Variables',nivar));
    catch
        nivar = ['mms' num2str(i) '_dis_numberdensity_brst'];
        nis = double(spdfcdfread(vfile,'Variables',nivar));
    end
    
    
    
    tv =spdfcdfread(vfile,'Variables','Epoch');
    tP =spdfcdfread(Pfile,'Variables','Epoch');
    tB = spdfcdfread(filename,'Variables','Epoch');
    tE = spdfcdfread(Efile,'Variables',tEvar);
    rGSE = zeros(size(BGSE));
    vGSE = rGSE; EGSE = rGSE;ueGSE=vGSE;
    EGSEf = smoo(EGSEf,100);
    for j = 1:3
        rGSE(:,j) = interp1(tr,rGSEslow(:,j),tB);
        vGSE(:,j) = interp1(tv,vGSEslow(:,j),tB);
        ueGSE(:,j) = interp1(tP,ueGSEslow(:,j),tB);
        EGSE(:,j) = interp1(tE(1:100:end),EGSEf(1:100:end,j),tB);
    end
    EGSE(:,4) = sqrt(EGSE(:,1).^2+EGSE(:,2).^2+EGSE(:,3).^2);
    rGSE(:,4) = sqrt(rGSE(:,1).^2+rGSE(:,2).^2+rGSE(:,3).^2);
    vGSE(:,4) = sqrt(vGSE(:,1).^2+vGSE(:,2).^2+vGSE(:,3).^2);
    ueGSE(:,4) = sqrt(ueGSE(:,1).^2+ueGSE(:,2).^2+ueGSE(:,3).^2);
    Pxx = interp1(tP,Pxxs,tB);
    Pxy = interp1(tP,Pxys,tB);
    Pxz = interp1(tP,Pxzs,tB);
    Pyy = interp1(tP,Pyys,tB);
    Pyz = interp1(tP,Pyzs,tB);
    Pzz = interp1(tP,Pzzs,tB);
    
    Tpar = interp1(tP,Tpars,tB);
    Tperp = interp1(tP,Tperps,tB);
    
    Tipar = interp1(tv,Tipars,tB);
    Tiperp = interp1(tv,Tiperps,tB);
    
    ne = interp1(tP,nes,tB);
    ni = interp1(tv,nis,tB);
    
    JGSE = 1e-6*([ni ni ni].*vGSE(:,1:3)-[ne ne ne].*ueGSE(:,1:3));
    
    %Rotate the pressure tensor to the magnetic field
    clear Ppar Pperp1 Pperp2 Pparp1 Pparp2 Pp1p2
    Ppar = zeros(size(Pxx));
    Pperp1 = Ppar;
    Pperp2 = Ppar;
    Pparp1 = Ppar;
    Pparp2 = Ppar;
    Pp1p2 = Ppar;
    
    % Normalizing the B-fields
    Bmod=BGSE(:,4);
    bxu=BGSE(:,1)./Bmod;
    byu=BGSE(:,2)./Bmod;
    bzu=BGSE(:,3)./Bmod;

    %make  unit vector perp to B in the xz-plane
    bp1x=bzu;
    bp1z=-bxu;
    Bmodp=sqrt(bp1x.^2+bp1z.^2);
    bp1x=bp1x./Bmodp;
    bp1z=bp1z./Bmodp;
    bp1y=bp1z*0;

    %make  unit vector perp to B and bp1
    bp2x=byu.*bp1z-bzu.*bp1y;
    bp2y=bzu.*bp1x-bxu.*bp1z;
    bp2z=bxu.*bp1y-byu.*bp1x;
    
    for tindex = 1:length(tB)
        R=[[bxu(tindex) byu(tindex)  bzu(tindex)] ;
            [bp1x(tindex) bp1y(tindex)  bp1z(tindex)] ;
            [bp2x(tindex) bp2y(tindex)  bp2z(tindex)] ] ;


        Pten=[[Pxx(tindex) Pxy(tindex) Pxz(tindex)];
            [Pxy(tindex) Pyy(tindex) Pyz(tindex)];
            [Pxz(tindex) Pyz(tindex) Pzz(tindex)]];

        Pten2=R*Pten*R';
        Ppar(tindex)=Pten2(1,1);
        A=Pten2(2,2);B=Pten2(3,3);C=Pten2(2,3);
        thetar = atan2(2*C,A-B)/2;
        R2 = [1,0,0;0,cos(thetar),sin(thetar);0,-sin(thetar),cos(thetar)];
        Ptenf = R2*Pten2*R2';
        Pperp1(tindex)= Ptenf(2,2);
        Pperp2(tindex)= Ptenf(3,3);
        Pparp1(tindex)=Ptenf(1,2);
        Pparp2(tindex)=Ptenf(1,3);
        Pp1p2(tindex)=Ptenf(2,3);
    end
    save(savefile,'EGSE','BGSE','tB','tr','tf','rGSE','vGSE','Pxx','Pxy','Pxz',...
        'Pyy','Pyz','Pzz','Ppar','Pperp1','Pperp2','Pparp1','Pparp2','Pp1p2',...
        'ne','ueGSE','ni','JGSE','Tpar','Tperp','Tipar','Tiperp');
    toc 
    tic
end
