function X=GeneratePeptideMedium(nth,nt)

aenv=4;
centre=1;

arad=2.250000;
len=51.060677;
Xa=GenerateRodSmoothCaps(nth,nt,arad,aenv,len,centre);
Xa=RotatePoints(Xa,[0 0 0],-pi/2,1);
Xa=RotatePointsAxisAngle(Xa,[0.793966 0.000000 -0.607962 0.520577]);
Xa=TranslatePoints(Xa,[2.401536 2.550737 -5.454571]);

arad=2.250000;
len=51.114220;
Xb=GenerateRodSmoothCaps(nth,nt,arad,aenv,len,centre);
Xb=RotatePoints(Xb,[0 0 0],-pi/2,1);
Xb=RotatePointsAxisAngle(Xb,[0.885352 0.000000 0.464921 0.582156]);
Xb=TranslatePoints(Xb,[3.646448 -1.305767 5.114625]);

arad=2.250000;
len=51.169350;
Xc=GenerateRodSmoothCaps(nth,nt,arad,aenv,len,centre);
Xc=RotatePoints(Xc,[0 0 0],-pi/2,1);
Xc=RotatePointsAxisAngle(Xc,[-0.002477 0.000000 0.999997 0.101977]);
Xc=TranslatePoints(Xc,[-6.670410 -0.490515 1.048739]);

X=MergeVectorGrids(Xa,Xb);
X=MergeVectorGrids(X,Xc);

