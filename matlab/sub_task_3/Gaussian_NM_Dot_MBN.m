m=[2 5 6 10:14 34:38 44:46 55:57 62];
r=[1:15 20:39 44:46 55:57 62];%44:46 49:57 62:71 74];
v0=permute(coordinates(46,r,:,37,37),[2 3 1 4 5]);
v0_norm=sqrt(trace(v0'*v0));
v0=v0/v0_norm;
vm0=permute(coordinates(46,m,:,37,37),[2 3 1 4 5]);
vm0=vm0/v0_norm;
Tr_v=zeros(72,72,100,'double');
Tr_vm=zeros(72,72,100,'double');
Tr_vm_v=zeros(72,72,100,'double');
Tr_v_candidates=zeros(72,72,100,'double');
Tr_vm_candidates=zeros(72,72,100,'double');
Tr_vm_v_candidates=zeros(72,72,100,'double');
for il=1:72
    for ir=1:72
        for iv=1:100
            v=permute(coordinates(iv,r,:,il,ir),[2 3 1 4 5]);
            v_norm=sqrt(trace(v'*v));
            v=v/v_norm;
            vm=permute(coordinates(iv,m,:,il,ir),[2 3 1 4 5]);
            vm=vm/v_norm;
            Tr_v(il,ir,iv)=abs(sqrt(trace(v'*v0)));
            Tr_vm(il,ir,iv)=abs(sqrt(trace(vm'*vm0)));
        end
        Tr_v_max=max(Tr_v(il,ir,:));
        Tr_v_candidates(il,ir,:)=Tr_v(il,ir,:)>=0.7*Tr_v_max;
        Tr_vm_max=max(Tr_vm(il,ir,:));
        Tr_vm_candidates(il,ir,:)=Tr_vm(il,ir,:)>0.7*Tr_vm_max;
    end
end




if false
    boopboop=Tr_vm./Tr_v;boopboop(boopboop>1)=0;boopboop=boopboop.*Tr_v_candidates;
    figure;
    for il=1:72
        subplot(2,2,1);contourf(abs(permute(Tr_v(il,:,:),[2 3 1])),20,'edgecolor','none');
        subplot(2,2,2);contourf((permute(frequencies(43:49,il,:),[3 1 2])),20,'edgecolor','none');
        subplot(2,2,3);contourf(abs(permute(Tr_vm(il,:,:),[2 3 1])),20,'edgecolor','none');
        subplot(2,2,4);contourf(abs(permute(boopboop(il,:,:).*Tr_v_candidates(il,:,:),[2 3 1])),20,'edgecolor','none');
        title(num2str(il));
        pause
    end
end

if false
    figure;subplot(3,2,1);contourf(abs((Tr_vm_candidates(33:41,33:41,44)+Tr_v_candidates(33:41,33:41,44)).*Tr_vm(33:41,33:41,44)),20,'edgecolor','none');
    subplot(3,2,2);contourf(abs((Tr_vm_candidates(33:41,33:41,45)+Tr_v_candidates(33:41,33:41,45)).*Tr_vm(33:41,33:41,45)),20,'edgecolor','none');
    subplot(3,2,3);contourf(abs((Tr_vm_candidates(33:41,33:41,46)+Tr_v_candidates(33:41,33:41,46)).*Tr_vm(33:41,33:41,46)),20,'edgecolor','none');
    subplot(3,2,5);contourf(abs((Tr_vm_candidates(33:41,33:41,47)+Tr_v_candidates(33:41,33:41,47)).*Tr_vm(33:41,33:41,47)),20,'edgecolor','none');
    subplot(3,2,6);contourf(abs((Tr_vm_candidates(33:41,33:41,48)+Tr_v_candidates(33:41,33:41,48)).*Tr_vm(33:41,33:41,48)),20,'edgecolor','none');
    figure;subplot(3,2,1);contourf((permute(frequencies(44,33:41,33:41),[2 3 1])),20,'edgecolor','none');colorbar;caxis([305 370]);
    subplot(3,2,2);contourf((permute(frequencies(45,33:41,33:41),[2 3 1])),20,'edgecolor','none');colorbar;caxis([305 370]);
    subplot(3,2,3);contourf((permute(frequencies(46,33:41,33:41),[2 3 1])),20,'edgecolor','none');colorbar;caxis([305 370]);
    subplot(3,2,5);contourf((permute(frequencies(47,33:41,33:41),[2 3 1])),20,'edgecolor','none');colorbar;caxis([305 370]);
    subplot(3,2,6);contourf((permute(frequencies(48,33:41,33:41),[2 3 1])),20,'edgecolor','none');colorbar;caxis([305 370]);
end

