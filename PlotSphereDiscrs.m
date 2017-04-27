function PlotSphereDiscrs(m,q,a,vis,wd,ht,ftnm,ftsz,filename)

% plots figure of Sphere Meshes with m and q subdivisions per face, a is
% radius

[x,X]=GenerateSpherePoints(m,q,a);
N=length(x)/3;Q=length(X)/3;

figure(1);clf;
[x1,x2,x3]=ExtractComponents(x);
plot3(x1,x2,x3,'c.');hold on;grid on;
if ~isempty(vis)
    [x1,x2,x3]=GenerateSphereSurfaceForVisualisation(vis(1),vis(2),a);
    hs=surf(x1,x2,x3,0*x3);shading flat;alpha(0.7);
    light;
end
box on;axis equal;
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

figure(1);clf;
[x1,x2,x3]=ExtractComponents(X);grid on;
plot3(x1,x2,x3,'m.');hold on;
if ~isempty(vis)
    [x1,x2,x3]=GenerateSphereSurfaceForVisualisation(vis(1),vis(2),a);
    hs=surf(x1,x2,x3,0*x3);shading flat;alpha(0.7);
    light;
end
box on;axis equal;
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
%print(1,'-dpdf',[filename num2str(Q) '.pdf']);
print(1,'-dpng','-r600',[filename num2str(Q) '.png']);

