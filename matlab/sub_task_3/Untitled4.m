DATA_46(DATA_46<0)=NaN;
temp=permute(DATA_46(1,30:44,27:42),[2 3 1]);
figure;contourf(L_prop(27:42),L_prop(30:44),temp);xlabel('\DeltaR');ylabel('\DeltaL');


figure;contourf(L_prop(:),L_prop(:),permute(DATA_46(1,:,:),[2 3 1]));xlabel('\DeltaR');ylabel('\DeltaL');
figure;plot(1:101,R_dihedrals-R_dihedrals(1),1:101,L_dihedrals-L_dihedrals(1))


v0=permute(coordinates(46,:,:,37,37),[2 3 1 4 5]);
v=permute(coordinates(:,:,:,37,37),[2 3 1 4 5]);
for il=1:72
    for ir=1:72
        for iv=1:222
            dotV(il,ir,iv)=sqrt(sum(sum(permute(coordinates(iv,:,:,il,ir),[2 3 1 4 5]).*v0)));
        end
    end
end
figure;contourf(abs(permute(dotV(:,37,:),[1 3 2])),'edgecolor','none')
I=[];
for il=1:72
    for ir=1:72
        [m,ii]=max(abs(dotV(il,ir,:)));
        ff(il,ir)=frequencies(ii,il,ir);
        I(il,ir)=ii;
    end
end
ff(ff<0)=NaN;
figure;mesh(L_prop(30:42),L_prop(30:44),ff(30:44,30:42));xlabel('\DeltaR');ylabel('\DeltaL');

for il=1:72
    for ir=1:72
        [m,ii]=max(real(dotV(il,ir,:)));
        ff(il,ir)=frequencies(ii,il,ir);
        I(il,ir)=ii;
    end
end
ff(ff<0)=NaN;
figure;mesh(L_prop(30:42),L_prop(30:44),ff(30:44,30:42));xlabel('\DeltaR');ylabel('\DeltaL');


