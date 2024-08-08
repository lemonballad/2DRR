% Get X Y Z from COMbH
XX=[(sum(X(resid==43))+sum(X(resid==45)))/(length(X(resid==43))+length(X(resid==45))) ...
(sum(X(resid==64))+sum(X(resid==68)))/(length(X(resid==64))+length(X(resid==68))) ...
(sum(X(resid==93))+sum(X(resid==97)))/(length(X(resid==93))+length(X(resid==93)))];
YY=[(sum(Y(resid==43))+sum(Y(resid==45)))/(length(Y(resid==43))+length(Y(resid==45))) ...
(sum(Y(resid==64))+sum(Y(resid==68)))/(length(Y(resid==64))+length(Y(resid==68))) ...
(sum(Y(resid==93))+sum(Y(resid==97)))/(length(Y(resid==93))+length(Y(resid==93)))];
ZZ=[(sum(Z(resid==43))+sum(Z(resid==45)))/(length(Z(resid==43))+length(Z(resid==45))) ...
(sum(Z(resid==64))+sum(Z(resid==68)))/(length(Z(resid==64))+length(Z(resid==68))) ...
(sum(Z(resid==93))+sum(Z(resid==97)))/(length(Z(resid==93))+length(Z(resid==93)))];

C=7*10^-3;
%XVG Extract for forces
Xff=C*ff(1,:)+XX;
Yff=C*ff(2,:)+YY;
Zff=C*ff(3,:)+ZZ;
C=7.5*10^-3;
%XVG Extract for forces
tipX=C*ff(1,:)+XX;
tipY=C*ff(2,:)+YY;
tipZ=C*ff(3,:)+ZZ;

[XX(1) YY(1) ZZ(1);Xff(1) Yff(1) Zff(1);tipX(1) tipY(1) tipZ(1)]
[XX(2) YY(2) ZZ(2);Xff(2) Yff(2) Zff(2);tipX(2) tipY(2) tipZ(2)]
[XX(3) YY(3) ZZ(3);Xff(3) Yff(3) Zff(3);tipX(3) tipY(3) tipZ(3)]