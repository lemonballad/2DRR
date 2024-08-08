range=1:2499;
t_axis=[];nan_axis=[];
for ii=0:19
    current_range=range+ii*2500;
    t_axis=[t_axis current_range];
    nan_axis=[nan_axis (ii+1)*2500];
end
t=1:50001;nan_axis=[nan_axis 50001];
test1=NaN(19122,50001);
test1(:,t_axis)=test;
mu_test1=mean(test,2);mu_test1=repmat(mu_test1,1,21);
test1(:,nan_axis)=mu_test1;

M=matfile('C:\Users\Thomas\Desktop\mbc\mbc_XYZ.mat');
test=M.cH2O_x;