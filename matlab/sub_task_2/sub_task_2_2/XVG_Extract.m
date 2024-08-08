r0=[3.094+3.000+3.063+3.150...
    3.832+3.568+3.557+3.820...
    2.824+2.816+2.545+2.553]/4;
r1=[3.043 3.705 2.677];
r2=[3.202 3.632 2.709];
r3=[3.313 3.632 2.697];
v1=[3.094-3.000 3.832-3.568 2.824-2.816];
v2=[3.094-3.000 3.832-3.557 2.824-2.545];
n=cross(v1,v2);n=n/norm(n);
dr=0.01;
r=-0.1:dr:0.1;r=r-0.036;
lr=length(r);


forces=zeros(lr,1586,3,'double');
for ii=1:lr
    xvg_file=['C:\Users\Thomas\Desktop\MBC_FOLDER\MBC_XVG\force_' num2str(ii) '.xvg'];
    forces(ii,:,1:3)=reshape(dlmread(xvg_file,'\t',4779,1),3,1586)';
end

hforce=forces(:,1538:end,:);
pforce=forces(:,1:1537,:);

hf0=permute(hforce(11,:,:),[2 3 1]);
hf0=permute(repmat(hf0,1,1,lr),[3 1 2]);
dhf=hforce-hf0;

pf0=permute(pforce(11,:,:),[2 3 1]);
pf0=permute(repmat(pf0,1,1,lr),[3 1 2]);
dpf=pforce-pf0;

pf_res=zeros(lr,length(fres)-1,3,'double');
dpf_res=zeros(lr,length(fres)-1,3,'double');
xp_res=zeros(lr,length(fres)-1,'double');
yp_res=zeros(lr,length(fres)-1,'double');
zp_res=zeros(lr,length(fres)-1,'double');
for ii=1:length(fres)-1;
    pf_res(:,ii,:)=sum(pforce(:,ires(ii):fres(ii),:),2);
    dpf_res(:,ii,:)=sum(dpf(:,ires(ii):fres(ii),:),2);
    xp_res(:,ii)=sum(xp(ires(ii):fres(ii)))/length(ires(ii):fres(ii));
    yp_res(:,ii)=sum(yp(ires(ii):fres(ii)))/length(ires(ii):fres(ii));
    zp_res(:,ii)=sum(zp(ires(ii):fres(ii)))/length(ires(ii):fres(ii));
end

labfont=16;lwidax=4;
figure;contourf(1:49,r,log10(sqrt(sum(dhf.^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;

figure;contourf(1:1537,r,log10(sqrt(sum(dpf.^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;

figure;contourf(1:49,r,log10(sqrt(sum(hforce.^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;

figure;contourf(1:1537,r,(sqrt(sum(pforce.^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;

figure;contourf(1:153,r,(sqrt(sum(pf_res.^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;

figure;contourf(1:153,r([1:10 12:21]),(sqrt(sum(dpf_res([1:10 12:21],:,:).^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;ylabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax,...
    'XTick',[43 68 93],'XTickLabel',res([43 68 93]));
title('\DeltaF (kJ mol^{-1} nm^{-1})');caxis([000 2000]);

figure;contourf(1:153,1:20,log10(sqrt(sum(pf_res([1:10 12:21],:,:).^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;ylabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax,...
    'XTick',[43 68 93],'XTickLabel',res([43 68 93]),...
    'YTick',[1 6 15 20],'YTickLabel',[-0.1 -0.05 0.05 0.1]);%r([1:10 12:21])
title('log_{10}|\DeltaF/F_0|');

figure;contourf(1:153,r([1:10 12:21]),(sqrt(sum(dpf_res([1:10 12:21],:,:).^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;ylabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax,...
    'XTick',[43 68 93],'XTickLabel',res([43 68 93]));
title('\DeltaF (kJ mol^{-1} nm^{-1})');

figure;contourf(1:153,1:20,log10(sqrt(sum(dpf_res([1:10 12:21],:,:).^2,3))),50,'edgecolor','none');
colormap jet;colorbar;axis square;ylabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax,...
    'XTick',[43 68 93],'XTickLabel',res([43 68 93]),...
    'YTick',[1 6 15 20],'YTickLabel',[-0.1 -0.05 0.05 0.1]);%r([1:10 12:21])
title('log_{10}|\DeltaF/F_0|');

figure;plot(r,log10(sqrt(sum(pf_res(:,[43 45 64 68 93 97],:).^2,3))),...
    'linewidth',4);
axis square;ylabel('Force (kJ mol^{-1} nm^{-1})');xlabel('\Deltar (nm)');
legend('PHE43','ARG45','HIS64','VAL68','HIS93','HIS97');
set(gca,'fontsize',labfont,'linewidth',lwidax);
% hold on
% plot(r,log10(sqrt((sum(hforce.^2,3)))),'linewidth',4);

pf40s=sum(pf_res(:,[43 45],:),2);
pf60s=sum(pf_res(:,[64 68],:),2);
pf90s=sum(pf_res(:,[93 97],:),2);
figure;plot(r,(sqrt(sum(pf40s.^2,3))),'b-',...
    r,(sqrt(sum(pf60s.^2,3))),'r-',...
    r,(sqrt(sum(pf90s.^2,3))),'g-',...
    'linewidth',4);
% hold on
% plot(r,(sqrt(sum(sum(hforce.^2,3),2))),'k-','linewidth',4);
axis square;ylabel('Force (kJ mol^{-1} nm^{-1})');xlabel('\Deltar (nm)');
% legend('PHE43-ARG45','HIS64-VAL68','HIS93-HIS97','HEME');
set(gca,'fontsize',labfont,'linewidth',lwidax);ylim([500 3000]);


figure;plot(r,(sqrt(sum(pf_res(:,:,:).^2,3))),'linewidth',4);
axis square;ylabel('Force (kJ mol^{-1} nm^{-1})');xlabel('\Deltar (nm)');
% legend('PHE43','ARG45','HIS64','VAL68','HIS93','HIS97');
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;scatter3(xp_res(1,[43 45 64 68 97]),...
    yp_res(1,[43 45 64 68 97]),...
    zp_res(1,[43 45 64 68 97]),100,'k.');
hold on;
for ii=1:4:lr
quiver3(xp_res(ii,[43 45 64 68 97]),...
    yp_res(ii,[43 45 64 68 97]),...
    zp_res(ii,[43 45 64 68 97]),...
    pf_res(ii,[43 45 64 68 97],1),...
    pf_res(ii,[43 45 64 68 97],2),...
    pf_res(ii,[43 45 64 68 97],3))
end
