po1=' SCF Done:  E(ROB3LYP) =  ';
p1=' SCF Done:  E(RB3LYP) =  ';
po2=' Input=MBO_';
p2=' Input=1MBN_';

[e1,d1]=rdEngGau('list1.txt',p1,p2);
[e2,d2]=rdEngGau('list2_2.txt',po1,po2);

E1=e1(e1~=0);
D1=d1(e1~=0);
E2=e2(e2~=0);
D2=d2(e2~=0);

if D1(1)>D1(end)
    D1=fliplr(D1);
    E1=fliplr(E1);
end

if D2(1)>D2(end)
    D2=fliplr(D2);
    E2=fliplr(E2);
end

E1=E1*219474.6313705;
E2=E2*219474.6313705;

figure
plot(D1,E1)
figure
plot(D2,E2)