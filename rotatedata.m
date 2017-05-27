function  fmat=rotatedata(thmat,phmat,Tht,Phi,Cntarr,Ne)

Tht=[Tht-180, Tht, Tht+180];
Phi=[Phi-360, Phi, Phi+360];

fmat=zeros([size(thmat) Ne]);

for l=1:Ne
   for k=1:size(thmat,1)
      Ctarr=squeeze(Cntarr(:,:,l));
      farr=[[Ctarr, Ctarr, Ctarr]; [Ctarr, Ctarr, Ctarr]; [Ctarr, Ctarr, Ctarr]];
      fmat(k,:,l)=interp2(Tht,Phi,farr,thmat(k,:),phmat(k,:));
   end
end