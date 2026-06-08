clearvars;

blockSize=0.2;

epsrel=0.01;
domain='i';

MColl=[3 6 12 24];
QColl=[12 24 48 96 192]; 

%MColl = [3 6 12];
%QColl = [12 24];

TablesForReport_Sphere_vs_cnearest(0.010,MColl,QColl);
