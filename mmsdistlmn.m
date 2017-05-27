function mmsdistlmn(directory)
cd(directory)
path(path,'C:\matlab_cdf361_patch');
path(path,directory);
Datadir = dir('*Data');
Datadir = [directory '\' Datadir.name];
load([Datadir '\MMS_lmn.mat']);
varlist = {'Bfl','Efl','gradnl','gradfl','perpplanegradl','perpplanegradln'};
for mms = 1:4
    distvars = dir([Datadir '\MMS' num2str(mms) '*ave*']);
    load([Datadir '\' distvars.name]);
    Bfl = Bf;
    perpplanegradn = squeeze(mean(gradfn(:,:,23:27,:),3));
    for tindex = 1:length(Bf(:,1))
        Bfl(tindex,1:3) = Bf(tindex,1:3)*Tlmn';
        Efl(tindex,:) = Ef(tindex,1:3)*Tlmn';
        gradnl(tindex,:) = gradn(tindex,:)*Tlmn';
        for enerind = 1:length(perpplanegrad(1,1,:))
            perpplanegradl(tindex,:,enerind)=squeeze(perpplanegrad(tindex,:,enerind))*Tlmn';
             perpplanegradln(tindex,:,enerind)=squeeze(perpplanegradn(tindex,:,enerind))*Tlmn';
            for thetaind = 1:length(gradf(1,1,:,1))
                gradfl(tindex,:,thetaind,enerind)=squeeze(gradf(tindex,:,thetaind,enerind))*Tlmn';
            end
        end
    end
    save([Datadir '\MMS' num2str(mms) 'ave'],varlist{:},'-append');
end
