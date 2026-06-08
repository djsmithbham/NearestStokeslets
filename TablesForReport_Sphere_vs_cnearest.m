function TablesForReport_Sphere_vs_cnearest(epsilon,MColl,QColl)

%Stokes' Law test case: this version with collpts and quadpts different and force points difference

%force discretised via nearest neighbour interpolation

% this version for direct comparison with cnearest/examples/example1sphere

U=[1,0,0];
Om=[1,0,0];
X0=[0,0,0];

a=1;
domain='i';
blockSize=0.2;

for iM=1:length(MColl)
    for iQ=1:length(QColl)
        tic;
        fprintf('iM = %i, iQ = %i\n',iM,iQ);
        [x,X]=GenerateSpherePoints(MColl(iM),QColl(iQ),a);
        [FF,MM,~] = SolveRigidResistanceMultipleModes(x,X,X0,[U' [0;0;0]],[[0;0;0] Om'],epsilon,domain,blockSize,'n'); 
        FTr = FF(:,1);
        MRo = MM(:,2);
        FTrError(iM,iQ)=(FTr(1)-6*pi*a)/6/pi/a;     % x1-components only for direct comparison with CNearest
        MRoError(iM,iQ)=(MRo(1)-8*pi*a^3)/8/pi/a^3;
        hf(iM)=CalcDiscr_h(x,blockSize);
        hq(iQ)=CalcDiscr_h(X,blockSize);
        Dof(iM)=length(x);
        NQ(iQ)=length(X)/3;
        NM(iM)=length(x)/3;
        walltime(iM,iQ)=toc;
    end
end

fid=fopen(['nearest_sphere_transrot_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'|    N     |    Q     |   Force Error   |  Moment Error   |    Time (s)     |\n');
fprintf(fid,'-----------------------------------------------------------------------------\n');
for iM=1:length(MColl)
    for iQ=1:length(QColl)
        fprintf(fid,['|' pad(num2str(NM(iM)),8,'left') '  |' pad(num2str(NQ(iQ)),8,'left') ...
                     '  |' pad(num2str(FTrError(iM,iQ),'%1.8f'),16,'left') ' |'  pad(num2str(MRoError(iM,iQ),'%1.8f'),16,'left') ...
                     ' |'  pad(num2str(walltime(iM,iQ),'%1.8f'),16,'left') ' |\n']);
    end
end
% 
% 
% 
%     fprintf(fid,'%i & %i &  ',NM(iM),Dof(iM));
%     fprintf(fid,'& %.3f',walltime(iM,:));
%     fprintf(fid,'\\\\ \n    ');
% end
% fprintf(fid,'\\end{tabular}\n');
fclose(fid);


fid=fopen(['translation_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Translating sphere, infinite fluid test\n');
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
    fprintf(fid,'& %1.4f',FTrError(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);

fid=fopen(['rotation_test_' num2str(epsilon) '.txt'],'w');
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
    fprintf(fid,'& %1.4f',MRoError(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);


