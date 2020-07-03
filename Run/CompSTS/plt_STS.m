clear; clc;
set(0, 'DefaultAxesFontSize',20);

[linestyles,MarkerEdgeColors,Markers,mycolor]=mystyles(20);
ccolor = distinguishable_colors(15);

FontAxis=24;
FontAxisLegend=26;
FontLegend=24;
FontText=20;
FontAnnotation=24;

FontAxisInset=22;
FontAxisLegendInset=22;
FontLegendInset=24;
FontTextInset=24;
FontAnnotationInset=20;
scrsz = get(groot,'ScreenSize');
ScreenTop=scrsz(4);
ScreenLeft=scrsz(3);

%close all
f=figure(1);
%clf

% Set the desired colormap.
%colormap jet;

% Set directories.
%Dir_data = '/Users/shukla/Research/Nice/STS_Tests/Run_256/FHpy/nrst200/CompSTS';
Dir_data = '../Run1/Run_STS';
Dir_plot = 'STS';
% File number.
nrst = 50;


%%======================================================
% Simulation parameters.

% Box length:
Lx = 2*pi; Ly = 2*pi;
%-------------------------------------
% Resolution:
Nx=256; Ny=256;
%-------------------------------------
% Spatial step size.
dx = Lx/Nx; dy = Ly/Ny;
%-------------------------------------
% Maximum wave number (after dealiasing).
kmax=floor(Nx/3);
%-------------------------------------
% No. of snapshots saved (k-slice).
Ntt = 4096;
%-------------------------------------
% Temporal step size.
dt = 2.0*10^(-4);
%-------------------------------------
% No. of steps between two consecutive snapshots (k-slice).
wt = 10;
%-------------------------------------
% In the forced dissipated simulation, mean density is not fixed to 1.
% Therefore, we have to find it first.
% Now as rho_0 is not 1, the xi and c in the GPE code have no meaning.
% The numbers alpha0 and g that appear in the equation are the only meaningful
% quantity, now they can be regaded as being aribitriry specified.
% xi and c must be computed in post facto.
% The only relation that hold is the c^2 xi^2 = 2 alpha0^2. 
% However, alpha0 and g must be computed using the xi and c from the code to begin with.
% inc_ == > Incorrect or as in the GPE code.
% Speed of sound.
csound_inc = 1;
%-------------------------------------
% Healing length.
nhl = 1;
xi_inc= nhl*dx;

% Assumed rho0 = 1.
rho0_inc = 1.0;%
alpha0 = csound_inc*xi_inc/sqrt(2);
g = alpha0/(xi_inc^2*rho0_inc);
%-------------------------------------
%-------------------------------------
% Compute actual rho0
a=load([Dir_data,'/data/norm.dat_nrst',num2str(nrst)]);
rho0 = mean(a(1:end,2))/(2*pi)^2;

% Speed of sound.
csound = sqrt(2*alpha0*g*rho0);
% Healing length
xi = sqrt(2)*alpha0/csound;

%b=load('../data/n0.dat_nrst22');
%b(:,2) = b(:,2)*(2*pi)^3
%-------------------------------------
%-------------------------------------
% Temporal length over which STS is computed.
T=dt*Ntt*wt;
%-------------------------------------
% Frequency step.
dw = 2*pi/T;
%-------------------------------------
% Maximum frequency from simulation.
wmaxS = (Ntt/2)*dw;
%-------------------------------------
%%======================================================

% Generate wave vector.
kx=[-Nx/2:1:Nx/2-1];
% Generate frequencies.
w=[-Ntt/2:1:Ntt/2-1]*dw;
% Set positive frequencies upto kmax.
k = [0:1:floor(Nx/3)];
%-------------------------------------

% Theoretical or fitted Dispersion relation.
wtheo = csound*abs(kx).*sqrt(1 + 0.5*xi^2*kx.^2);
%wtheo = csound*abs(kx).*sqrt(0 + 0.5*xi^2*kx.^2);
%-------------------------------------
% Maximum frequency from theory.
wmax = csound*abs(kmax).*sqrt(1 + 0.5*xi^2*kmax.^2);

%%======================================================

R = load([Dir_data,'/CompSTS/','sts2d_nrst',num2str(nrst),'.dat']);

% Rearrange in the wave number space. (like fftshift)
A(Ntt/2+1:Ntt,:)=R(1:Ntt/2,:);
for i = Ntt/2+1:Ntt;
A(i-Ntt/2,:) = R(i,:);
end;


% Normalize at each wave number.
B = max(A,[],1);
C = A./B;

% Plots
h(1) = pcolor(xi*k,(xi/csound)*w,C);shading interp;colorbar;caxis([0 1]);
hold on;
h(2) = plot(xi*kx,(xi/csound)*(-wtheo - 4.1),'Color','w','linestyle','--','linewidth',1);
%h(3) = plot(xi*kx,0.5*(xi/csound)*(-wtheo - 4.1),'Color','w','linestyle','--','linewidth',1);
%h(2) = plot(kx,-wtheo-23.45,'Color','w','linestyle','--','linewidth',1);
%plot(kx,wtheo-23.45,'Color','w','linewidth',1);xlim([0 floor(Nx/3)]);
%hold off;
%ylim([-wmax wmax]);

xlim([0 2.1]);
ylim([-4 2]);

%% Zoom
%xlim([0 3]);
%ylim([-6 2]);


%text(70,-10,'$$\omega(k)\sim k^2$$','Interpreter','Latex','Color','w','Fontsize',FontText);
xlabel('$\xi k$','FontSize',FontAxisLegend,'Interpreter','Latex');
ylabel('$\xi \omega(k)/c$','FontSize',FontAxisLegend,'Interpreter','Latex');

set(gca,'xtick',[0 0.5 1.0 1.5 2.0]);
set(gca,'xTickLabel',{'0';'0.5';'1.0';'1.5';'2.0'});
set(gca,'ytick',[-4 -2 0 2]);
set(gca,'yTickLabel',{'-4';'-2';'0';'2'});


Leg=cell(1);
Leg{1}=[''];
Leg{2}=['$\omega_{B}$'];
%Leg{3}=['$\frac{1}{2}\omega_{B}$'];

legend([h],Leg,'FontSize',FontLegend,'Interpreter','Latex','Textcolor','w','Box','off');

text(0.25,-3,'$\xi k_{\rm max}=2.09$','Interpreter','Latex','Color','w','FontSize',FontText);

%print ('-depsc2',[Dir_plot,'/','STS_R1']);
%print ('-dpdf',[Dir_plot,'/','STS_R1']);
print('-dpng',[Dir_plot,'/','STS_R1']);
savefig([Dir_plot,'/','STS_R1']);
