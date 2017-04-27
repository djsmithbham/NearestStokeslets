% Refinement heuristic applied to rotation of a unit sphere

clear all
close all

U=[1,0,0];
Om=[1,0,0];
X0=[0,0,0];

a=1;
domain='i';
blockSize=0.2;

epsilon=0.01;

MColl=[3 6];
QColl=[48 96];

for iM=1:length(MColl)
    for iQ=1:length(QColl)
        fprintf('iM = %i, iQ = %i',iM,iQ);
        [x,X]=GenerateSpherePoints(MColl(iM),QColl(iQ),a);
        [FRo,MRo,fRo,~]=SolveRigidResistance(x,X,X0,[0,0,0],Om,    epsilon,domain,blockSize);
        MRoStore(iM,iQ)=MRo(1);
        hf(iM)=CalcDiscr_h(x,blockSize);
        hq(iQ)=CalcDiscr_h(X,blockSize);
        Dof(iM)=length(x);
        NQ(iQ)=length(X)/3;
        NM(iM)=length(x)/3;
    end
end

fid=fopen(['rotation_test_alg1_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Rotating sphere, infinite fluid test\n');
fprintf(fid,'%%Radius, %f\n',a);
fprintf(fid,'%%Regularisation parameter, %f\n',epsilon);
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
    fprintf(fid,'& %1.4f',MRoStore(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)
