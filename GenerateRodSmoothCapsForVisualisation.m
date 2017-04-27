function [X1,X2,X3]=GenerateRodSmoothCapsForVisualisation(nth,nt,arad,aenv,len,centre)

%Envelope
E=@(t)(t<len-aenv).*(t>aenv)+(t<=aenv).*sqrt(1-((aenv-t)/aenv).^2)+(t>=len-aenv).*sqrt(1-((aenv-(len-t))/aenv).^2);

th=linspace(0,2*pi,nth+1)';

tt=linspace(0,len,nt+2);

X1=arad*cos(th)*E(tt);
X2=arad*sin(th)*E(tt);
X3=ones(nth+1,1)*tt;

switch centre
    case 0
        X3=X3;
    case 1
        X3=X3-len*0.5;
end

