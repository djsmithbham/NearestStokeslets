function TablesForReport_Prolate(ep,MColl,QColl)

%prolate spheroid test case

%force discretised via nearest neighbour interpolation

U=[1,0,0];
Om=[1,0,0];
X0=[0,0,0];

a=5.0;
c=1.0;
e=sqrt(a^2-c^2)/a;

d=[1;0;0];

L=log((1+e)/(1-e));

XA=8/3*e^3/(-2*e+(1+e^2)*L);
YA=16/3*e^3/(2*e+(3*e^2-1)*L);
XC=4/3*e^3*(1-e^2)/(2*e-(1-e^2)*L);
YC=4/3*e^3*(2-e^2)/(-2*e+(1+e^2)*L);

A_an=6*pi*a*(XA*d*d'+YA*(eye(3)-d*d'));
B_an=zeros(3);
C_an=8*pi*a^3*(XC*d*d'+YC*(eye(3)-d*d'));

R_an=[A_an B_an'; B_an C_an];

domain='i';
blockSize=0.2;

for iM=1:length(MColl)
    for iQ=1:length(QColl)
        tic;
        fprintf('iM = %i, iQ = %i\n',iM,iQ);
        [x,X]=GenerateSpherePoints(MColl(iM),QColl(iQ),1);
        [x1,x2,x3]=ExtractComponents(x);
        [X1,X2,X3]=ExtractComponents(X);
        x1=a*x1; X1=a*X1;
        x2=c*x2; X2=c*X2;
        x3=c*x3; X3=c*X3;
        x=[x1;x2;x3]; X=[X1;X2;X3];
        A  = CalcTransResistanceMatrix    (x,X,   ep,domain,blockSize);
        B  = CalcTransMomentCouplingMatrix(x,X,X0,ep,domain,blockSize);
        BT = CalcRotForceCouplingMatrix   (x,X,X0,ep,domain,blockSize);
        C  = CalcRotResistanceMatrix      (x,X,X0,ep,domain,blockSize);
        R=[A BT; B C];
        RError(iM,iQ)=norm(R-R_an)/norm(R_an);
        hf(iM)=CalcDiscr_h(x,blockSize);
        hq(iQ)=CalcDiscr_h(X,blockSize);
        Dof(iM)=length(x);
        NQ(iQ)=length(X)/3;
        NM(iM)=length(x)/3;
    end
end

fid=fopen(['prolate_test_' num2str(ep) '.txt'],'w');
fprintf(fid,'%%Resistance matrix for prolate spheroid, infinite fluid test\n');
fprintf(fid,'%%Long semi-axis, %f\n',a);
fprintf(fid,'%%Short semi-axis, %f\n',c);
fprintf(fid,'%%Regularisation parameter, %f\n',ep);
fprintf(fid,'    \\begin{tabular}{cccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &           & \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'          &     &           & \\(h_q \\)');fprintf(fid,'& %1.4f',hq);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF & \\(h_f\\) &          \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           &          \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i & %1.4f & ',NM(iM),Dof(iM),hf(iM));
    fprintf(fid,'& %1.4f',RError(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)

