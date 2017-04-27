function [Xa1,Xa2,Xa3,Xb1,Xb2,Xb3,Xc1,Xc2,Xc3]=GeneratePeptideSurfaceForVisualisation(nth,nt)

% Mesh functions are used only for visualisation, not computation

aenv=4;
centre=1;

arad=2.250000;
len=51.060677;
[Xa1,Xa2,Xa3]=GenerateRodSmoothCapsForVisualisation(nth,nt,arad,aenv,len,centre);
[Xa1,Xa2,Xa3]=RotateMesh(Xa1,Xa2,Xa3,[0 0 0],-pi/2,1);
[Xa1,Xa2,Xa3]=RotateMeshAxisAngle(Xa1,Xa2,Xa3,[0.793966 0.000000 -0.607962 0.520577]);
[Xa1,Xa2,Xa3]=TranslateMesh(Xa1,Xa2,Xa3,[2.401536 2.550737 -5.454571]);

arad=2.250000;
len=51.114220;
[Xb1,Xb2,Xb3]=GenerateRodSmoothCapsForVisualisation(nth,nt,arad,aenv,len,centre);
[Xb1,Xb2,Xb3]=RotateMesh(Xb1,Xb2,Xb3,[0 0 0],-pi/2,1);
[Xb1,Xb2,Xb3]=RotateMeshAxisAngle(Xb1,Xb2,Xb3,[0.885352 0.000000 0.464921 0.582156]);
[Xb1,Xb2,Xb3]=TranslateMesh(Xb1,Xb2,Xb3,[3.646448 -1.305767 5.114625]);

arad=2.250000;
len=51.169350;
[Xc1,Xc2,Xc3]=GenerateRodSmoothCapsForVisualisation(nth,nt,arad,aenv,len,centre);
[Xc1,Xc2,Xc3]=RotateMesh(Xc1,Xc2,Xc3,[0 0 0],-pi/2,1);
[Xc1,Xc2,Xc3]=RotateMeshAxisAngle(Xc1,Xc2,Xc3,[-0.002477 0.000000 0.999997 0.101977]);
[Xc1,Xc2,Xc3]=TranslateMesh(Xc1,Xc2,Xc3,[-6.670410 -0.490515 1.048739]);
