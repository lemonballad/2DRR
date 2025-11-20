m=[10:12 31:33];%[10:12 35:37];
% m=[2 5 6 10:13 15 34:37 39];%[10:12 35:37]; 1:15 16:19 20:39 40:43
r=[1:15 20:39];%44:46 49:57 62:71 74];
v0=permute(coordinates(46,r,:,37,37),[2 3 1 4 5]);
v0_norm=sqrt(trace(v0'*v0));
v0=v0/v0_norm;
% vm0=v0(m,:);
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
            v=v/1/v_norm;
%             vm=v(m,:);
            vm=permute(coordinates(iv,m,:,il,ir),[2 3 1 4 5]);
            vm=vm/v_norm;
            Tr_v(il,ir,iv)=abs(sqrt(trace(v'*v0)));
            Tr_vm(il,ir,iv)=abs(sqrt(trace(vm'*vm0)));
%             Tr_vm_v(il,ir,iv)=abs(sqrt(trace(vm'*vm)./trace(v'*v)));
        end
        Tr_v_max=max(Tr_v(il,ir,:));
        Tr_v_candidates(il,ir,:)=Tr_v(il,ir,:)>=0.8*Tr_v_max;
        Tr_vm_max=max(Tr_vm(il,ir,:));
        Tr_vm_candidates(il,ir,:)=Tr_vm(il,ir,:)>0.8*Tr_vm_max;
%         Tr_vm_v_max=max(Tr_vm_v(il,ir,:));
%         Tr_vm_v_candidates(il,ir,:)=Tr_vm_v(il,ir,:)>0.7*Tr_vm_v_max;
    end
end
