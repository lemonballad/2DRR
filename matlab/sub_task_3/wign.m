c=0.00003;hbar=5300;
m=1;
w=10;
x0=-50:50;
p0=x0/10^5/3;
n=1;
t=0:100;
for ix=1:101
    for ip=1:101
        for it=1:101
            x(ix,ip,it)=m*2*pi*c*w*x0(ix)*cos(2*pi*c*w*t(it))-p0(ip)*sin(2*pi*c*w*t(it));
            p(ix,ip,it)=p0(ip)*cos(2*pi*c*w*t(it))+2*pi*c*w*m*x0(ix)*sin(2*pi*c*w*t(it));
        end
    end
end
u=m*(2*pi*c*w)^2*x.^2/2+p.^2/2/m;
F=(-1)^n/pi/hbar*laguerreL(n,4*u/hbar/w).*exp(-2*u/hbar/w);
contourf(F(:,:,1),25,'edgecolor','none');

for it=1:101
    contourf(F(:,:,it),25,'edgecolor','none');
    pause
end