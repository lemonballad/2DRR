po1=' SCF Done:  E(UB3LYP) =  ';
p1=' SCF Done:  E(RB3LYP) =  ';
po2=' Input=1MBO_';
p2=' Input=1MBN_';

[e1,d1]=rdEngGau('list.txt',p1,p2);
[e2,d2]=rdEngGau('list2.txt',po1,po2);

E1=e1(e1~=0);
D1=d1(e1~=0);
E2=e2(e2~=0);
D2=d2(e2~=0);

figure
plot(D1,E1)
figure
plot(D2,E2)