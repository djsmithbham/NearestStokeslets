% sphere

clear all;

m=4;
q=10;
a=1;
vis=[100,50];
wd=7;
ht=6;
ftnm='times';
ftsz=10;
filename='sphereDiscrFig';

%%

PlotSphereDiscrs(m,q,a,vis,wd,ht,ftnm,ftsz,filename);

%%

% Peptide structure

% was 9,61 / 19,201 - need to systematically check

m=[6,21];
q=[11,51];

filename='peptideDiscrFig';

PlotPeptideDiscrs(m,q,vis,wd,ht,ftnm,ftsz,filename);
