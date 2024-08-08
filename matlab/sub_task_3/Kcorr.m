
n=length(mbcSS);indices=1:n/1;
traj=mbdSS(indices)-mean(mbdSS);temp=xcorr(traj);
ntraj=length(traj);
t_corr=[temp(ntraj:end)/ntraj/4 zeros(1,n-ntraj,'double')];

k_vc=(5308)^-2*trapz(t_corr)*4;figure;plot(t,t_corr)
title([num2str(k_vc) ' fs^{-1} ' num2str(1/k_vc) ' fs'])

