function [X1,X2,X3]=TranslateMesh(X1,X2,X3,a)

% performs translation on mesh X1,X2,X3 by 3-vector a for visualisation purposes

X1=X1+a(1);
X2=X2+a(2);
X3=X3+a(3);

