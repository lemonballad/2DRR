dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_map.txt',mbo(1).map(122:222,137:237));
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_phi_L.txt',(-50:50)-94);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_phi_R.txt',(-35:65)-109);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_freq.txt',mbo(2).freq(1:2501));
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_time.txt',mbo(2).time(1:2501));
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_omega_axis.txt',mbo(1).omega);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_spec_den.txt',mbo(1).spec_den);

dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_map.txt',mbn(1).map(122:222,122:222));
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_phi_L.txt',(-50:50)-81);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_phi_R.txt',(-50:50)-81);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_freq.txt',mbn(1).freq(42000:44500));
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_time.txt',(0:2500)*0.02);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_omega_axis.txt',mbn(1).omega);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_spec_den.txt',mbn(1).spec_den);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_LDA.txt',mbn(1).LDA-81);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbn_RDA.txt',mbn(1).RDA-81);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_LDA.txt',mbo(1).LDA-94);
dlmwrite('C:\Users\Thomas\Desktop\Myo_Fig\mbo_RDA.txt',mbo(1).RDA-109);







test=dlmread('C:\Users\Thomas\Desktop\Myo_Fig\mbo.txt');
testL=dlmread('C:\Users\Thomas\Desktop\Myo_Fig\phi_L.txt');
testR=dlmread('C:\Users\Thomas\Desktop\Myo_Fig\phi_R.txt');
figure;contourf(testL,testR,test,50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\Phi_R');ylabel('\Phi_L');


figure;plot(mbn(1).omega,mbn(1).spec_den,mbn(2).omega,mbn(2).spec_den,'Linewidth',3);
legend('filtered data','fit exponentials');xlim([0 100]);xlabel('\omega (cm^-^1)');ylabel('<\omega(0)|\omega(t)>');

plot(mbn(1).RDA);


plot((0:2501)*0.02,mbn(1).freq(42000:44501));xlim([0 50]);

figure;contourf((-35:65)-109,(-50:50)-94,test(122:222,137:237),50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\Phi_R');ylabel('\Phi_L');

figure;contourf((-50:50),(-50:50),mbn(1).map(122:222,122:222),50,'edgecolor','none');
colormap('jet');colorbar;xlabel('\DeltaR');ylabel('\DeltaL');