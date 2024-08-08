function y=expnfitt(x,t)
%reads the input parameters and generates curve g(t)
h=x';
if length(h)==3
for z=1:length(t)
     g(1,z)=(h(1)*exp(-t(z,1)/h(2)))+h(3);
end
    y=g;
end

if length(h)==5
for z=1:length(t)
     g(1,z)=h(1)*exp(-t(z,1)/h(2))+h(3)*exp(-t(z,1)/h(4))+h(5);
end
    y=g;
end

if length(h)==7
for z=1:length(t)
     g(1,z)=h(1)*exp(-t(z,1)/h(2))+h(3)*exp(-t(z,1)/h(4))+h(5)*exp(-t(z,1)/h(6))+h(7);
end
    y=g;
end

return


