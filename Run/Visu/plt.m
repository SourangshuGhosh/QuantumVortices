clear;clc;clf;
colormap Jet;
set(0,'DefaultAxesFontSize',15);
N = 128;

for t = 101:1:300;

fid=fopen(['rho/rho',num2str(t),'.dat'],'r');
dum = fread(fid,1,'float32');
rho = fread(fid,[N N],'float64');
dum = fread(fid,1,'float32');
fclose(fid);

rho = transpose(rho);

x=0:2*pi/N:2*pi-2*pi/N;
y=0:2*pi/N:2*pi-2*pi/N;

pcolor(x,y,rho);shading interp;colorbar;
caxis([0. 1.2]);
xlabel('x','fontsize',14,'fontweight','b')
ylabel('y','fontsize',14,'fontweight','b')
title(['\rho, Time = ',num2str(t*01)]);
axis('square');



%set(gca,'xtick',[])
%set(gca,'ytick',[])


print('-dpng',['plots/t',num2str(t,'%0.5d')]);
 
%clear;
clc;
clf;
end;
%clear k;
