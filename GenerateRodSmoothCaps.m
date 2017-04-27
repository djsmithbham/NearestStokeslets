function X=GenerateRodSmoothCaps(nth,nt,arad,aenv,len,centre)

%Envelope
E=@(t)(t<len-aenv).*(t>aenv)+(t<=aenv).*sqrt(1-((aenv-t)/aenv).^2)+(t>=len-aenv).*sqrt(1-((aenv-(len-t))/aenv).^2);

th=linspace(0,2*pi,nth+1)';
th(end)=[];

tt=linspace(0,len,nt+2);
tt(end)=[];
tt(1)=[];

X1=arad*cos(th)*E(tt);
X2=arad*sin(th)*E(tt);
X3=ones(nth,1)*tt;

X1=X1(:);
X2=X2(:);

% end points
X1=[X1; 0; 0  ];
X2=[X2; 0; 0  ];

switch centre
    case 0
        X3=X3(:);
        X3=[X3; 0; len];
    case 1
        X3=X3(:)-len*0.5;
        X3=[X3; -len*0.5; len*0.5];
end

X=[X1;X2;X3];
