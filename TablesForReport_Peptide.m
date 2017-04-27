function TablesForReport_Peptide(hCoarse,HCoarse,mlev,qlev,lenscale,epsrel,domain,blockSize)

[arad,aenv,len]=ApproxPeptideRodParams;

hTarget=len/2*2.^(-(0:mlev-1))*hCoarse; % [hCoarse hCoarse/2 hCoarse/4];
HTarget=len/2*2.^(-(0:qlev-1))*HCoarse; % [HCoarse HCoarse/2 HCoarse/4];

epsilon=epsrel*len/2;

kB=1.38064852e-23;  % Boltzmann constant in SI units
T=310;              % absolute temperature in Kelvin
mu=1e-3;            % dynamic viscosity in Pa.s

for iM=1:length(hTarget)
    for iQ=1:length(HTarget)
        tic;
        fprintf('iM = %i, iQ = %i',iM,iQ);
        [nth,nt,h]=CalcRodDiscretisationParameters(hTarget(iM),arad,aenv,len,blockSize);
        [Nth,Nt,H]=CalcRodDiscretisationParameters(HTarget(iQ),arad,aenv,len,blockSize);
        x=GeneratePeptideMedium(nth,nt);
        X=GeneratePeptideMedium(Nth,Nt);
        D=CalcFreeDiffusionTensor(x,X,[0 0 0],kB,mu,T,lenscale,epsilon,domain,blockSize);
        Dr=D(4:6,4:6);
        if sum(isnan(Dr))==0
            lambda1=min(eig(Dr));
        else
            lambda1=NaN;
        end
        tauStore(iM,iQ)=1/6/lambda1;
        hf(iM)=CalcDiscr_h(x,blockSize);
        hq(iQ)=CalcDiscr_h(X,blockSize);
        Dof(iM)=length(x);
        NQ(iQ)=length(X)/3;
        NM(iM)=length(x)/3;
        walltime(iM,iQ)=toc;
    end
end

%%
fid=fopen(['peptide_rotational_timescale_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Rotational diffusion timescale of peptide, infinite fluid test, results in nanoseconds\n');
fprintf(fid,'%%Regularisation parameter, %f\n',epsilon);
fprintf(fid,'    \\begin{tabular}{cccc|');
for iQ=1:length(HTarget)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &           & \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'          &     &           & \\(h_q \\)');fprintf(fid,'& %1.4f',hq);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF & \\(h_f\\) &          \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           &          \\\\[-0.8em] \n    ');
for iM=1:length(hTarget)
    fprintf(fid,'%i & %i & %1.4f & ',NM(iM),Dof(iM),hf(iM));
    fprintf(fid,'& %1.4f',tauStore(iM,:)*1e9);
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);

fid=fopen(['walltime_peptide_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'    \\begin{tabular}{ccc|');
for iQ=1:length(HTarget)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &    \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF &           \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           \\\\[-0.8em] \n    ');
for iM=1:length(hTarget)
    fprintf(fid,'%i & %i &  ',NM(iM),Dof(iM));
    fprintf(fid,'& %.3f',walltime(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);
