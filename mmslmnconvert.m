function mmslmnconvert(directory,startsec,endsec,varargin)
tic
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];

if nargin>3
    if isa(varargin{1},'double')
        l = varargin{1};
        m = varargin{2};
        n = varargin{3};
        l = l/norm(l);
        m = m/norm(m);
        n = n/norm(n);
        sz = size(l);
        if sz(1) == 1
            l = l';
            m = m';
            n = n';
        end
        Tlmn = [l';m';n'];
    else
        ERROR %Add new cases later
    end
else
    datestring = dir('mms1_fgm*');
    datestring = datestring.name;
    datestring = strsplit(datestring,'_');
    datestring = datestring{5};

    yr = str2double(datestring(1:4));
    mn = str2double(datestring(5:6));
    dy = str2double(datestring(7:8));
    h0 = str2double(datestring(9:10));
    m0 = str2double(datestring(11:12));
    s0 = str2double(datestring(13:14));

    sstart = mod(s0 + startsec,60);
    send = mod(s0 + endsec,60);
    mstart = mod(m0+floor((s0+startsec)/60),60);
    mend = mod(m0+floor((s0+endsec)/60),60);
    hrstart = mod(h0+floor((m0+floor((s0+startsec)/60))/60),24);
    hrend = mod(h0+floor((m0+floor((s0+endsec)/60))/60),24);

    B = [];
    v = [];
    for i = 1:4
        loadfile = [Datadir '\MMS' num2str(i) '_GSE'];
        load(loadfile);
        B = [B;tB,rGSE(:,1:3),BGSE];
        v = [v;tB,tB,vGSE];
    end


    [rlmn,~]=GSE_mva(B,v,yr,mn,dy,hrstart,mstart,sstart,hrend,mend,send,'lmax');
    l = rlmn(:,1);
    m = rlmn(:,2);
    n = rlmn(:,3);

    Tlmn = [l';m';n'];
end
savefile = [Datadir '\MMS_lmn'];
save(savefile, 'l','m','n','Tlmn');
for i = 1:4
    loadfile = [Datadir '\MMS' num2str(i) '_GSE'];
    load(loadfile);
    Blmn = zeros(size(BGSE));
    rlmn = zeros(length(BGSE(:,1)),3); vlmn = rlmn; Elmn = rlmn;
    for tindex = 1:length(BGSE(:,1))
        Blmn(tindex,1:3) = BGSE(tindex,1:3)*Tlmn';
        rlmn(tindex,:) = rGSE(tindex,1:3)*Tlmn';
        vlmn(tindex,:) = vGSE(tindex,1:3)*Tlmn';
        uelmn(tindex,:) = ueGSE(tindex,1:3)*Tlmn';
        Elmn(tindex,:) = EGSE(tindex,1:3)*Tlmn';
        Jlmn(tindex,:) = JGSE(tindex,1:3)*Tlmn';
    end
    Blmn(:,4) = BGSE(:,4);
    savefile = [Datadir '\MMS' num2str(i) '_lmn'];
    save(savefile, 'tB','Blmn','vlmn','rlmn','Elmn','ne','Ppar','Pperp1','Pperp2','uelmn','Jlmn','Tpar','Tperp')
end
loadfile = [Datadir '\MMS_FD'];
load(loadfile);
Gradnl = zeros(size(Gradn));
Gradperpnl = zeros(size(Gradperpn));
for tindex = 1:length(Gradn(:,1))
    Gradnl(tindex,:) = Gradn(tindex,:)*Tlmn';
    Gradperpnl(tindex,:) = Gradperpn(tindex,:)*Tlmn';
end
save(loadfile, 'Gradnl', 'Gradperpnl','-append')