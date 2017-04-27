function [nth,nt,h]=CalcRodDiscretisationParameters(hTarget,arad,aenv,len,blockSize)

% to assist with discretisation of rods...
% for given target max distance to nearest neighbour hTarget, this function
% estimates the required discretisation parameters nth, nt

nth=ceil(2*pi*arad/hTarget);
nt=ceil(len/hTarget);
centre=0;
x=GenerateRodSmoothCaps(nth,nt,arad,aenv,len,centre);
[h] = CalcDiscr_h(x,blockSize);