function PlotPeptideDiscrs(m,q,vis,wd,ht,ftnm,ftsz,filename);

x=GeneratePeptideMedium(m(1),m(2));
X=GeneratePeptideMedium(q(1),q(2));
N=length(x)/3;Q=length(X)/3;

figure(1);clf;
[x1,x2,x3]=ExtractComponents(x);
if ~isempty(vis)
    [Xa1,Xa2,Xa3,Xb1,Xb2,Xb3,Xc1,Xc2,Xc3]=GeneratePeptideSurfaceForVisualisation(vis(1),vis(2));
    hs=surf(Xa1,Xa2,Xa3,0*Xa3);shading flat;alpha(0.7);hold on;
    hs=surf(Xb1,Xb2,Xb3,0*Xb3);shading flat;alpha(0.7);
    hs=surf(Xc1,Xc2,Xc3,0*Xc3);shading flat;alpha(0.7);
    light;
end
plot3(x1,x2,x3,'c.');view([-37.5 30]);axis([-12 12 -25 25 -20 20]);
box on;axis equal;grid on;
set(gca,'tickdir','out');
hx=xlabel('\(x_1\)','interpreter','latex');
hy=ylabel('\(x_2\)','interpreter','latex');
hz=zlabel('\(x_3\)','interpreter','latex');
set(hx,'fontsize',ftsz);set(hx,'fontname',ftnm);
set(hy,'fontsize',ftsz);set(hy,'fontname',ftnm);
set(hz,'fontsize',ftsz);set(hz,'fontname',ftnm);
set(gca,'fontsize',ftsz);set(gca,'fontname',ftnm);
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
%print(1,'-dpdf',[filename num2str(N) '.pdf']);
print(1,'-dpng','-r600',[filename num2str(N) '.png']);

%
figure(1);clf;hold on;
[x1,x2,x3]=ExtractComponents(X);
if ~isempty(vis)
    [Xa1,Xa2,Xa3,Xb1,Xb2,Xb3,Xc1,Xc2,Xc3]=GeneratePeptideSurfaceForVisualisation(vis(1),vis(2));
    hs=surf(Xa1,Xa2,Xa3,0*Xa3);shading flat;alpha(0.7);hold on;
    hs=surf(Xb1,Xb2,Xb3,0*Xb3);shading flat;alpha(0.7);
    hs=surf(Xc1,Xc2,Xc3,0*Xc3);shading flat;alpha(0.7);
    light;
end
plot3(x1,x2,x3,'m.');view([-37.5 30]);axis([-12 12 -25 25 -20 20]);
box on;axis equal;grid on;
set(gca,'tickdir','out');
hx=xlabel('\(x_1\)','interpreter','latex');
hy=ylabel('\(x_2\)','interpreter','latex');
hz=zlabel('\(x_3\)','interpreter','latex');
set(hx,'fontsize',ftsz);set(hx,'fontname',ftnm);
set(hy,'fontsize',ftsz);set(hy,'fontname',ftnm);
set(hz,'fontsize',ftsz);set(hz,'fontname',ftnm);
set(gca,'fontsize',ftsz);set(gca,'fontname',ftnm);
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
%print(1,'-dpdf',[filename num2str(N) '.pdf']);
print(1,'-dpng','-r600',[filename num2str(Q) '.png']);

