%% Read files
if false
% HEME
mHEME=matfile('C:\Users\Thomas\Desktop\mbc\HEME\mbc_XYZ.mat','Writable',true);
m.cHEM_x=[];m.cHEM_y=[];m.cHEM_z=[];
carts_str='ATOM\s+\d+\s+\S+\s+HEM\s+155\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';

parfor ifile=0:100
    M=matfile(['C:\Users\Thomas\Desktop\mbc\HEME\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    M.cHEM_x=[];M.cHEM_y=[];M.cHEM_z=[];
    filePDB=['C:\Users\Thomas\Desktop\mbc\HEME\mbc_HEMEs' num2str(ifile) '.pdb'];
    f_in=fopen(filePDB);
    line=fgetl(f_in);
    cHEM_x=zeros(49,3000,'double');%19122
    cHEM_y=zeros(49,3000,'double');
    cHEM_z=zeros(49,3000,'double');
    iatom=0;
    itime=1;
    while ischar(line)
        if regexp(line,'ATOM');
            iatom=iatom+1;
            carts_data=regexp(line,carts_str,'names');
            cHEM_x(iatom,itime)=str2double(carts_data.x);
            cHEM_y(iatom,itime)=str2double(carts_data.y);
            cHEM_z(iatom,itime)=str2double(carts_data.z);
            if iatom==49
                iatom=0;
                itime=itime+1;
            end;
            if mod(itime,100)==0 && iatom==1;
                [ifile itime]
            end;
        end
        line=fgetl(f_in);
    end
    fclose(f_in);
    cHEM_x(:,itime:end)=[];
    cHEM_y(:,itime:end)=[];
    cHEM_z(:,itime:end)=[];
    M.cHEM_x=[M.cHEM_x cHEM_x];
    M.cHEM_y=[M.cHEM_y cHEM_y];
    M.cHEM_z=[M.cHEM_z cHEM_z];
end
for ifile=0:100
    M=matfile(['C:\Users\Thomas\Desktop\mbc\HEME\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    m.cHEM_x=[m.cHEM_x M.cHEM_x];
    m.cHEM_y=[m.cHEM_y M.cHEM_y];
    m.cHEM_z=[m.cHEM_z M.cHEM_z];
    ifile
end

% PROT
mPROT=matfile('C:\Users\Thomas\Desktop\mbc\PROT\mbc_XYZ.mat','Writable',true);
m.cPROT_x=[];m.cPROT_y=[];m.cPROT_z=[];
carts_str='ATOM\s+\d+\s+\S+\s+\S+\s+\d+\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';
parfor ifile=0:100
    M=matfile(['C:\Users\Thomas\Desktop\mbc\PROT\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    M.cPROT_x=[];M.cPROT_y=[];M.cPROT_z=[];
    filePDB=['C:\Users\Thomas\Desktop\mbc\PROT\mbc_PROTs' num2str(ifile) '.pdb'];
    f_in=fopen(filePDB);
    line=fgetl(f_in);
    cPROT_x=zeros(1537,3000,'double');% mbc: 1537 mbd: 
    cPROT_y=zeros(1537,3000,'double');
    cPROT_z=zeros(1537,3000,'double');
    iatom=0;
    itime=1;
    while ischar(line)
        if regexp(line,'ATOM');
            iatom=iatom+1;
            carts_data=regexp(line,carts_str,'names');
            cPROT_x(iatom,itime)=str2double(carts_data.x);
            cPROT_y(iatom,itime)=str2double(carts_data.y);
            cPROT_z(iatom,itime)=str2double(carts_data.z);
            if iatom==1537
                iatom=0;
                itime=itime+1;
            end;
            if mod(itime,100)==0 && iatom==1;
                [ifile itime]
            end;
        end
        line=fgetl(f_in);
    end
    fclose(f_in);
    cPROT_x(:,itime:end)=[];
    cPROT_y(:,itime:end)=[];
    cPROT_z(:,itime:end)=[];
    M.cPROT_x=[M.cPROT_x cPROT_x];
    M.cPROT_y=[M.cPROT_y cPROT_y];
    M.cPROT_z=[M.cPROT_z cPROT_z];
end
for ifile=0:100
    M=matfile(['C:\Users\Thomas\Desktop\mbc\PROT\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    m.cPROT_x=[m.cPROT_x M.cPROT_x];
    m.cPROT_y=[m.cPROT_y M.cPROT_y];
    m.cPROT_z=[m.cPROT_z M.cPROT_z];
    ifile
end

% H2O
m=matfile('C:\Users\Thomas\Desktop\mbc\H2O\mbc_XYZ.mat','Writable',true);
m.cH2O_x=[];m.cH2O_y=[];m.cH2O_z=[];
carts_str='ATOM\s+\d+\s+\S+\s+\S+\s+\d+\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';
parfor ifile=60:69
    M=matfile(['C:\Users\Thomas\Desktop\mbc\H2O\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    M.cH2O_x=[];M.cH2O_y=[];M.cH2O_z=[];
    filePDB=['C:\Users\Thomas\Desktop\mbc\H2O\mbc_H2Os' num2str(ifile) '.pdb'];
    f_in=fopen(filePDB);
    line=fgetl(f_in);
    cH2O_x=zeros(19122,3000,'double');% mbd: 19119 mbc: 19122
    cH2O_y=zeros(19122,3000,'double');
    cH2O_z=zeros(19122,3000,'double');
    iatom=0;
    itime=1;
    while ischar(line)
        if regexp(line,'ATOM');
            iatom=iatom+1;
            carts_data=regexp(line,carts_str,'names');
            cH2O_x(iatom,itime)=str2double(carts_data.x);
            cH2O_y(iatom,itime)=str2double(carts_data.y);
            cH2O_z(iatom,itime)=str2double(carts_data.z);
            if iatom==19122
                iatom=0;
                itime=itime+1;
            end;
            if mod(itime,100)==0 && iatom==1;
                [ifile itime]
            end;
        end
        line=fgetl(f_in);
    end
    fclose(f_in);
    cH2O_x(:,itime:end)=[];
    cH2O_y(:,itime:end)=[];
    cH2O_z(:,itime:end)=[];
    M.cH2O_x=[M.cH2O_x cH2O_x];
    M.cH2O_y=[M.cH2O_y cH2O_y];
    M.cH2O_z=[M.cH2O_z cH2O_z];
end
for ifile=30:69
    M=matfile(['C:\Users\Thomas\Desktop\mbc\H2O\mbc_XYZ' num2str(ifile) '.mat'],'Writable',true);
    m.cH2O_x=[m.cH2O_x M.cH2O_x];
    m.cH2O_y=[m.cH2O_y M.cH2O_y];
    m.cH2O_z=[m.cH2O_z M.cH2O_z];
    ifile
end
end

%% Deal with Charges
if false %% remove atoms with 0 charge  
    cHEM_q=repmat(mHEME.q,1,74970);
    cPROT_q=repmat(mPROT.q,1,74970);
    cH2O_q=repmat(mH2O.q,1,74970);

    cPROT_x(cPROT_q(:,1)==0,:)=[];
    cPROT_y(cPROT_q(:,1)==0,:)=[];
    cPROT_z(cPROT_q(:,1)==0,:)=[];
    cPROT_q(cPROT_q(:,1)==0,:)=[];
    
    cHEM_x(cHEM_q(:,1)==0,:)=[];
    cHEM_y(cHEM_q(:,1)==0,:)=[];
    cHEM_z(cHEM_q(:,1)==0,:)=[];
    cHEM_q(cHEM_q(:,1)==0,:)=[];
    
end

%% Calculate Solute-Solvent interaction
if false
nt=size(cHEM_x,2);
Solu_Solv=zeros(1,nt,'double');
epsilon_0=8.85*10^-12; % In SI units
for iatom=1:size(cHEM_x,1)
Solu_Solv=Solu_Solv+sum((1/4/pi/epsilon_0)*repmat(cHEM_q(iatom,:),size(cH2O_q,1),1).*cH2O_q./...
        sqrt(...
        (cH2O_x-repmat(cHEM_x(iatom,:),size(cH2O_x,1),1)).^2+...
        (cH2O_y-repmat(cHEM_y(iatom,:),size(cH2O_y,1),1)).^2+...
        (cH2O_z-repmat(cHEM_z(iatom,:),size(cH2O_z,1),1)).^2),1);

Solu_Solv=Solu_Solv+sum((1/4/pi/epsilon_0)*repmat(cHEM_q(iatom,:),size(cPROT_q,1),1).*cPROT_q./...
        sqrt(...
        (cPROT_x-repmat(cHEM_x(iatom,:),size(cPROT_x,1),1)).^2+...
        (cPROT_y-repmat(cHEM_y(iatom,:),size(cPROT_y,1),1)).^2+...
        (cPROT_z-repmat(cHEM_z(iatom,:),size(cPROT_z,1),1)).^2),1);
    iatom
end

Solu_Solv=Solu_Solv*10^9; % Accounting for the denominator needing to go from nanometers to meters
Solu_Solv=Solu_Solv*(1.6*10^-19)^2; % Converting the charge from e- to coulombs
Solu_Solv=Solu_Solv*5.034*10^22; % Convert from joules to wavenumbers

t=0:0.004:(size(cHEM_x,2)-1)*0.004;
figure;plot(t,Solu_Solv);
end
% M=matfile('C:\Users\Thomas\Desktop\mbc\mbc_XYZ.mat','Writable',true);Solu_Solv_mbc=M.Solu_Solv;
% M=matfile('C:\Users\Thomas\Desktop\mbd\mbd_XYZ.mat','Writable',true);Solu_Solv_mbd=M.Solu_Solv;

%% Plot
if false
    for ii=0:1
    labfont=16;lwidax=2;
    n=1;
    trange=1:n:125000;%mbc:1:25001    mbd: 1:25001
    if 0
        molecule='COMb';
        color=[0 0 1];
        traj=mbcSS;
    else
        molecule='DeoxyMb';
        color=[1 0 0];
        traj=mbdSS;
    end
    Solu_Solv_trj=traj(trange);
    dt=4*n;
    t=0:dt/1000/n:(length(traj)-1)*dt/1000/n;
    Nf=length(Solu_Solv_trj);
    omega=((-1/2:1/Nf:1/2-1/Nf)+1/Nf/2)*2*pi/dt/0.00003;
    mean_int=mean(Solu_Solv_trj);
    var_int=var(Solu_Solv_trj);
    std_int=sqrt(var_int);
    [mean_int std_int]
    figure;plot(t(trange),Solu_Solv_trj-mean_int,'Color',color);xlabel('t (ps)');ylabel('Solute-Solvent Energy (cm^{-1})');title(molecule);
    set(gca,'fontsize',labfont,'linewidth',lwidax);axis square;xlim([min(t(trange)) max(t(trange))]);
    temp=xcorr(Solu_Solv_trj-mean_int);%,'coeff');
    t_corr=temp(Nf:end);
    time=t(trange)-t(trange(1));
    f=fit(time',t_corr','exp2','StartPoint',[0.2 -4 0.8 -12])%,'Lower',[-0.076 -3.748],'Upper',[-0.071 -0.7]);
%     f=fit(time',t_corr','exp1','StartPoint',[0.2 -4]);%,'Lower',[-0.076 -3.748],'Upper',[-0.071 -0.7]);
%     apod=f.a*apodfun(time,0,0,-f.b,-f.b,false)+f.c*apodfun(time,0,0,-f.d,-f.d,false);figure;plot(time,t_corr,time,t_corr.*apod,time,apod);
    apod=apodfun(time,0,0,300,2,false);figure;plot(time,t_corr,time,t_corr.*apod,time,apod);%mbc: ~15 mbd: ~11
    f_corr=abs(fftshift(fft(t_corr.*apod)));f_corr=f_corr-f_corr(end);
    figure;semilogx(1000*time,t_corr,'Color',color,'linewidth',lwidax);xlabel('t (fs)');ylabel('C(t)');title(molecule);
    set(gca,'fontsize',labfont,'linewidth',lwidax,'XTick',[10 100 1000 10^4]);axis square;xlim(1000*[min(t(trange)) max(t(trange))]);ylim([-inf 1]);
    
    figure;plot(time,t_corr,'Color',color,'linewidth',lwidax);xlabel('t (ps)');ylabel('C(t)');title(molecule);
    set(gca,'fontsize',labfont,'linewidth',lwidax);axis square;xlim([min(t(trange)) max(t(trange))]);ylim([-inf 1]);
    
    figure;semilogx(omega,f_corr/max(f_corr),'linewidth',lwidax,'Color',color);
    xlim([1 500]);ylim([0 inf]);xlabel('\omega (cm^{-1})');ylabel('C(\omega)');title(molecule);
    set(gca,'fontsize',labfont,'linewidth',lwidax,'XTick',[1 10 100 1000]);axis square;
    
    figure;plot(omega,f_corr/max(f_corr),'linewidth',lwidax,'Color',color);
    xlim([0 250]);ylim([0 inf]);xlabel('\omega (cm^{-1})');ylabel('C(\omega)');title(molecule);
    set(gca,'fontsize',labfont,'linewidth',lwidax);axis square;
    end
 %     annotation('textbox',[0.35 0.7 0.45 0.2],'String',{['\mu_0 = ' num2str(mean_int,'%.4f') ' cm^{-1}'];['\Delta = \pm' num2str(std_int,'%.2f') ' cm^{-1}']},...
%     'FontSize',labfont,'fontweight','normal','backgroundcolor','white','color','black');
   
%     cent_mass=sum(omega(Nf/2+1:end).*(f_corr(Nf/2+1:end).^2))/sum((f_corr(Nf/2+1:end).^2))
    cent_mass=sum(omega(Nf/2+1:end).*abs(f_corr(Nf/2+1:end)))/sum(abs(f_corr(Nf/2+1:end)))
end

