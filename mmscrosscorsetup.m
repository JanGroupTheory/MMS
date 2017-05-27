function mmscrosscorsetup(directory)
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];
for i = 1:4
        loadfile = [Datadir '\MMS' num2str(i) 'ave.mat'];
        load(loadfile);
        Fzp = squeeze(mean(Fave(1:5,:,:),1));
        Fzm = squeeze(mean(Fave(46:50,:,:),1));
        Fp = squeeze(mean(Fave(23:27,:,:),1));
        frontname = ['Fzp' num2str(i)];
        perpname = ['Fp' num2str(i)];
        backname = ['Fzm' num2str(i)];
        eval([backname '= Fzm;']);
        eval([perpname '= Fp;']);
        eval([frontname '= Fzp;']);
        eval(['B' num2str(i) ' = Bfl;']);
end
varlist = {'Fzm1','Fzm2','Fzm3','Fzm4','Fzp1','Fzp2','Fzp3','Fzp4','Fp1','Fp2','Fp3','Fp4','B1','B2','B3','B4'};
save([Datadir '\MMSxcor.mat'],varlist{:});