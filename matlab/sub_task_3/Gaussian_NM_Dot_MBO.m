m=[1 12:17 36:41 46:48 61:63];
r=[1 4:17 22:41 46:48 61:63];
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
    bleepbleep=Tr_vm./Tr_v;bleepbleep(bleepbleep>1)=0;
    figure;
    for il=1:72
        subplot(2,2,1);contourf(abs(permute(Tr_v(il,:,:),[2 3 1])),20,'edgecolor','none');
        subplot(2,2,2);contourf((permute(frequencies(43:49,il,:),[3 1 2])),20,'edgecolor','none');
        subplot(2,2,3);contourf(abs(permute(Tr_vm(il,:,:),[2 3 1])),20,'edgecolor','none');
        subplot(2,2,4);contourf(abs(permute(bleepbleep(il,:,:).*Tr_v_candidates(il,:,:),[2 3 1])),20,'edgecolor','none');colorbar;
        title(num2str(il));
        pause
    end
end
