clear all;
close all;

%% --------------------------------------------------------------

% All figures

FiguresForReport;

TablesForReport_Heuristic;

%% --------------------------------------------------------------

% Peptide tables

% generate discretisations in angstrom (length is 50 angstrom)
lenscale=1e-10;


blockSize=0.2;

epsrel=0.01;
domain='i';

hCoarse=0.125;       % coarsest discretisation for force
HCoarse=hCoarse/4;   % coarsest discretisation for quadrature

mlev=3;
qlev=3;

TablesForReport_Peptide(hCoarse,HCoarse,mlev,qlev,lenscale,epsrel,domain,blockSize);

%% --------------------------------------------------------------

% Sphere tables

TablesForReport_Heuristic;

MColl=[3 6 12 24];
QColl=[12 24 48 96 192]; 

TablesForReport_Sphere(0.020,MColl,QColl);
TablesForReport_Sphere(0.010,MColl,QColl);
TablesForReport_Sphere(0.005,MColl,QColl);

%%
TablesForReport_Sphere(1e-6,MColl,QColl);
