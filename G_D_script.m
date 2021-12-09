clear all 

global Imax Jmax x xo y
global alphaF palat mp0 m2 Z 
global gmax P0 reminN reminS D0
global JR betaL betaN betaSi betaF
global L N0 S0 dil
global toppanelpos botpanelpos
toppanelpos = [.13,.51,.78,.4];
botpanelpos = [.13,.09,.78,.4];

%%%%%PARAMETERS%%%%%
%Carbon,Nitrogen,Si density
cb = 0.11e-6; % [µgC/µm^3] carbon content of cytoplasm (from Sthratmann,1967)
cm = 0.6e-6;  % [µgC/µm^3] carbon content of membranes
gb = 7;       % C:N ratio of cytoplasm [µgC/µgN]
gm = 38;      % C:N ratio of membrane [µgC/µgN]
Sishell = 10.3E-6; % silicate concentration in diatom shell [µg Si µm-3]
%cL = 0.0318*(gb/cb)^(-2/3);   % light affinity coefficient [µgC/day/(W m^-2)/µm^2]

%Thickness
h = 0.008;    % cell membranes [µm]
t = 0.01;     % shell [µm] 

%Affinity coefficients
kappa = 0.2;     % light affinity investment return [µm^-1] (also called absorption coefficient)
cN = 2.89E-4*(gb/cb)^(-1/3); % nutrient affinity coefficient [L/day/µm]
cP = 0.227*(gb)^(-1); % food affinity coefficient [L/day/µgC]
eP = 0.3408*(gb)^(-1) ; % nutrient affinity investment return [L/day/µgC]
alphaN = 1.25E-6;  % L day-1 (µg C)-1 µm2 (0.68) or L d-1 µm2
alphaL = 0.18 * 1E-11 * 24 * 3600*1.875;   % quantum yield (scaled)[µgC day-1 (µE m-2 s-1)-1 µm-2] 
qyield = 0.07; % quantum yield {µg C / µmol photons PAR]

%Food 
alphaF = 2E-8; %0.018 (L d-1 ?)
palat = 0.6; % diatom palatibility factor
epsilon = 7; % predator-prey size ratio (linear dimension)
width = 1; % in orders of magnitude
mp0 = 0.1; % grazing constant 
m2 = 0.01; % quadratic viral lysis [(µC L^-1)^-1 d^-1]
gmax = 1; %maximal growth rate [d^-1]
Z0 = 0.06; %copepods [µgC L^-1] (predation by higher trophic levels)
reminS = 0.8; reminN = 0.1;

%Metabolic cost 
JR = 0.04; % basal respiration [µgC µgC^-1 day^-1] (beta0)
betaL = 0.08;   % cost of light harvesting  [µgC µgC^-1]
betaN = 0.45;   % cost of nitrogen uptake  [µgC µgN^-1]
betaSi = 0.45;  % cost of silicate uptake  [µgC µgSi^-1]
betaF = 0.45;   % cost of particulate matter uptake  [µgC µgC^-1]

%Classes
Imax = 21; %Number of size classes 
Jmax = 12; %Number of classes for vacuolation or trophic strategies 

%Biomass 
Binit = 1;  % Initial total biomass abundance [mgC m^-3] or µgC m^-3 ??? 
BGtoD = 0.7;  % biomass split BGtoD in generalists, 1- BGtoD in diatoms

rstar = 2;  % crossover from r3 tp r absolute nutrient affinity 
xstar = 4*pi/3*rstar^3*cb; % mass equivalence for rstar

%%%%%VARIABLES%%%%%
%===primary trait space: x -> cell size (carbon content) and v -> vacuole ===
xo = 10.^linspace(-7,1,Imax);
[x,y] = meshgrid(xo,linspace(0,1,Jmax)); 
vo = linspace(0,1,Jmax+1); % vacuole size for diatoms -> C to membranes and 
[x,v] = meshgrid(xo,vo(1:Jmax));
r = (3/(4*pi*cb)*x./(1-v)).^(1/3);
%Investment in biomass photosynthesis
phib = 0; 
%Carbon amount apportioned to membrane mass
xm = 4*pi*r.^2*h*cm.*(1 + v.^(2/3));
phim = xm./x; 

%Generalists
G.x = x; G.r = ones(Jmax,1)*r(1,:); %(only radius with v=0)
G.phim = ones(Jmax,1)*phim(1,:); %phi_M(x,0) (first line of phim corresponds to v=0) 
phil = 1 - phim - phib - y; phil(phil<0) = 0; 
phif = y; phif(phil<0) = 0;
G.phil = phil;
G.phif = phif;

%Diatoms
D.x=x;D.v=v;D.r=r;
D.phim=phim;
D.phim(D.phim > 1) = 1;
phil = 1 - phim - phib; phil(phil<0) = 0;
D.phil = phil;

D.shell = 2.28E-3*(D.r).^(0.72); % thickness of diatom shell (t=1.62V^0.24)
D.s = 4*pi*Sishell*(D.r.^2).*D.shell; % diatom Si requirement [µg Si]
D.n = xm/gm + (D.x - D.phim.*D.x)/gb; % diatom N requirement [µg N]
%D.n(D.n<0)=0;
G.n = G.phim.*G.x/gm + (G.x-G.phim.*G.x)/gb; 

% C:Si ratio
D.CS = D.x./D.s; 
%C:N ratio 
D.CN= D.x./D.n;
G.CN = G.x./G.n; 

%===Affinities===
D.AL = alphaL*pi*(D.r.^2).*(1 - exp(-kappa*D.phil.*D.r.*(1-D.v)))./D.x; 
D.AN = 4/3*pi*cb*alphaN*(D.r)./(1 + (D.x/xstar).^(-2/3))./D.x; 
D.AN(D.phim == 1) = 0;


G.AL = alphaL*pi*(G.r.^2).*(1 - exp(-kappa*G.phil.*G.r))./G.x;
G.AN = alphaN*((G.r).^(-2))./(1 + (G.r/rstar).^(-2));
G.AF = alphaF*G.phif;


G.JF = zeros(size(G.AF));

%===Trophic interactions===
Theta = @(rpredator,rprey) exp(-(log(epsilon*rprey./rpredator)).^2/(2*width^2)); %represents size preference for smaller cells (theta_i,j)
%creates a function handle to an anonymous function 

T = zeros(Jmax,Imax,Imax);  
for i = 1:Imax  % calculated once 
    T(:,:,i) = Theta(r(1,i),D.r); %3D matrix with all possible predator prey interactions (prey are diatoms ?) 
end

%%%%%DYNAMICS%%%%%
%===Time===
tf = 365;
tspan=[1:tf];
years = 2;
tspan_y = [1:tf*years];

%===Seasonal or steady state===
seasonal = 1;

%===Initial environmental conditions===
if seasonal==1
    L = 105*(1-cos(2*pi*tspan/365))+25; L=L'; 
    dil = 0.9*(1+(-L+min(L))/210)+0.1;
%       dil= ones(tf,1);
%       dil_1 = 0.45*(1+cos(2*pi*tspan'/90))+0.1;
%       dil(89:134)=dil_1(1:46);
%       dil(135:200)= 0.1;
%       dil(201:246)=dil_1(45:90);

    Z=ones(tf,1)*5; Z(1:75)=5; Z(76:135)=0.25*tspan(76:135)-13.75; Z(136:239)=5/104*tspan(136:239)+13.5; Z(240:269)=-0.6*tspan(240:269)+168.4; Z(270:365)=-2/51*tspan(270:365)+5+2*320/51;   
    %Z = Z*1000 ;
%       Z_1 = 2.5*(1+cos(2*pi*tspan/60))+5;
%       Z_2= 1.25*(1+cos(2*pi*tspan/60))+5;
%       Z(61:121)=Z_1(30:90);
%       Z(214:274)=Z_2(30:90);
    
    L_y =[];dil_y=[];Z_y=[];
    for i=1:years
        L_y=[L_y;L]; dil_y = [dil_y;dil]; Z_y = [Z_y;Z];
    end 
    L=L_y;dil = dil_y; Z=Z_y;
else %steady state 
    dil = 0.5; %dil
    L=50; % light L [µmol photons m^-2 s^-1]
    Z=5; 
end 
    
N0 = 1200;  % nitrate concentration
S0 = 800; % silicate concentration Si [mgN m^-3]

P0 = ones([Jmax Imax])*0.001*N0; %Background P concentration [µgC L^-1] (to change)
D0 = P0;
for i=1:Imax*Jmax
    if D.AN(i)==0
        D0(i)=0;
    end 
end

D.b = ones(size(phil))*Binit*(1-BGtoD)./(Imax*Jmax); %initial total biomass abundance for diatoms 
G.b = ones(size(phil))*Binit*BGtoD./(Imax*Jmax); %initial total biomass abundance for generalists
D.b0 = reshape(D.b,[Jmax*Imax 1]); %initial biomass, nutrients at the end
G.b0 = reshape(G.b,[Jmax*Imax 1]);
B0 = [D.b0;G.b0];
maxind = 2*Imax*Jmax; 
if seasonal==1
    B0 = [D.b0;G.b0;N0;S0]; %2*Jmax*Imax + 2 vector
    maxind = Jmax*2*Imax+2;
end 

%===Solving the equation===
options = odeset('NonNegative',[1:maxind],'RelTol',1e-4,'AbsTol',1e-4); 
[t C] = ode45(@(t,C) equdyn(t,C,D,G,r,T,seasonal),tspan_y,B0,options);
[No_use HD HG] = equdyn(t,C,D,G,r,T,seasonal)  ; %output of D and G
D = HD; G = HG;

D.B = zeros(tf*years,Jmax,Imax); 
G.B = zeros(tf*years,Jmax,Imax);
D.BB = D.B;
G.BB = D.BB;
D.Psi = D.BB;
G.Psi = D.BB; 
D.Simp = zeros(tf*years,Jmax);
D.Shan = D.Simp ;
G.Simp = D.Simp; G.Shan = D.Shan;
for it=1:length(t)
    D.B(it,:,:) = reshape(squeeze(C(it,1:Jmax*Imax)),[Jmax Imax]); 
    G.B(it,:,:) = reshape(squeeze(C(it,Jmax*Imax+1:Jmax*2*Imax)),[Jmax Imax]);
    hd = reshape(D.B(it,:,:),[Jmax Imax]); %   for diversity index calculation
    hg = reshape(G.B(it,:,:),[Jmax Imax]);
    D.BB(it,:,:) = sum(hd(:));
    G.BB(it,:,:) = sum(hg(:));
    D.Psi(it,:,:) = D.B(it,:,:)./D.BB(it,:,:);
    G.Psi(it,:,:) = G.B(it,:,:)./G.BB(it,:,:);
    for j=1:Jmax
        D.Simp(it,j) = sum(D.Psi(it,j,:));
        G.Simp(it,j) =  sum(G.Psi(it,j,:));
        D.Shan(it,j) = - sum(D.Psi(it,j,:).*log(D.Psi(it,j,:)));
        G.Shan(it,j) = - sum(G.Psi(it,j,:).*log(G.Psi(it,j,:)));
    end 
    
end

if seasonal==1
     N = C(:,2*Jmax*Imax+1);
     Si = C(:,2*Jmax*Imax+2); 
end 
 
toppanelpos = [.13,.51,.78,.4];
botpanelpos = [.13,.09,.78,.4];

%Plot steady state 
if seasonal==0 
    figure(1)
    clf;
    filename = 'teststeadygeneralists.gif';
    for tt=1:10:tf
         futG=reshape(G.B(tt,:,:),[Jmax Imax]);
         imagesc(log10(xo),1-y(:,1),futG);
         set(gca, 'YDir', 'normal');
         colormap(jet(50));
         c=colorbar;
         ylabel(c, 'P \rm [mgC L^{-1} ]','FontSize',16)
         ylabel('\rm phi_L','FontSize',16);title('Evolution of total biomass for generalists','FontSize',16);
         subtitle(sprintf('t=%f days',tt-1));
         xlabel('\it x \rm [µg C]','FontSize',16);
         drawnow
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
         if tt == 1;
               imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
         else
               imwrite(imind,cm,filename,'gif','WriteMode','append');
         end
     end 

    figure(2)
    clf;
    filename = 'teststeadydiatoms2.gif';
    for tt=1:10:tf
        futD=reshape(D.B(tt,:,:),[Jmax Imax]);
        imagesc(log10(xo),D.v(:,1),futD);
        set(gca, 'YDir', 'normal');
        colormap(jet(50));
        c=colorbar;
        ylabel(c, 'P \rm [mgC L^{-1} ]','FontSize',16)
        ylabel('\it v','FontSize',16);title('Evolution of total biomass for diatoms','FontSize',16);
        subtitle(sprintf('t=%f days',tt-1));
        xlabel('\it x \rm [µg C]','FontSize',16);
        drawnow
        frame = getframe(2);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if tt == 1;
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end 

    figure(3)
    D.Btot = squeeze(sum(D.B,2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.B,2));
    plot(log10(xo),D.Btot(end,:),log10(xo),G.Btot(end,:));
     legend('diatoms','generalists');
    xlabel('\it x \rm [µg C]','FontSize',24);
    ylabel('Biomass \rm [mgC.m^{-3}]');
else 
 
    % Plot seasonal cycle 
    D.Btot = squeeze(sum(D.B,2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.B,2));


    figure(1)
    %padded=padarray(squeeze(D.Btot(:,:).*1e-3).',[1,1],'post');
    %pcolor(D.Btot);
    D.Btot_y=D.Btot((years-1)*tf+1:years*tf,:);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),D.Btot_y'); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf])
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    colormap(jet(50));
    C=colorbar;
    title('Biomass evolution for diatoms ');
    xlabel('Months')
    xlabel(C,'Biomass \rm [µgC.m^{-3}]');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp   
    ij=0;
    for size_lab=1:2:Imax
        ij=ij+1;
        size_labels(ij,:)=num2str(xo(size_lab),'%.1e');
    end
    yticklabels(size_labels)

    figure(2)
    %padded=padarray(squeeze(D.Btot(:,:).*1e-3).',[1,1],'post');
    %pcolor(D.Btot);
    G.Btot_y=G.Btot((years-1)*tf+1:years*tf,:);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),G.Btot_y'); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf])
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    colormap(jet(50));
    C=colorbar;
    title('Biomass evolution for generalists ');
    xlabel('Months')
    xlabel(C,'Biomass \rm [µgC.m^{-3}]');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp   
    ij=0;
    for size_lab=1:2:Imax
        ij=ij+1;
        size_labels(ij,:)=num2str(xo(size_lab),'%.1e');
    end
    yticklabels(size_labels)

    %nutrients
    figure(3)
    plot(t((years-1)*tf+1:years*tf),N((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),Si((years-1)*tf+1:years*tf))
    legend('N \rm [\mugN m^{-3}]','Si \rm [\mugS m^{-3}]');
    title('Evolution of nutrient concentration'); 
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')

    D.Btrait = zeros(tf*years,Imax);
    for it=1:tf*years
        for i=1:Imax
            k=max(D.B(it,:,i));
            kbis= max(find(D.B(it,:,i)==k));
            D.Btrait(it,i) = D.v(kbis,1);
        end 
    end 

    G.Btrait = zeros(tf*years,Imax);
    for it=1:tf*years
        for i=1:Imax
            k=max(G.B(it,:,i));
            kbis= max(find(G.B(it,:,i)==k));
            G.Btrait(it,i) = 1-y(kbis,1); 
        end 
    end 


    figure(4)
    D.Btrait_y=D.Btrait((years-1)*tf+1:years*tf,:);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),D.Btrait_y'); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf]);
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    C=colorbar;
    title('Dominant trait values for diatoms ');
    xlabel('Months')
    xlabel(C,'\it v (%)');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp 
    % 
    figure(5)
    G.Btrait_y=G.Btrait((years-1)*tf+1:years*tf,:);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),G.Btrait_y'); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf]);
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    C=colorbar;
    title('Dominant trait values for generalists ');
    xlabel('Months')
    xlabel(C,'\it \phi_L');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp 
    
    figure(6)
    D.Btotot = sum(D.Btot,2);
    G.Btotot = sum(G.Btot,2);
    plot(t((years-1)*tf+1:years*tf),D.Btotot((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),G.Btot((years-1)*tf+1:years*tf)); 
    legend('diatoms','generalists');
    ylabel('Biomass [mgC.m^{-3}]');
    title('Seasonal variation of plankton biomass'); 
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')
    
    figure(7)
    clf;
    subplot('Position',toppanelpos)
    i1 = imagesc(t((years-1)*tf+1:years*tf),D.v(:,1),D.Simp((years-1)*tf+1:years*tf)'); set(gca, 'YDir', 'normal');
    colormap(jet(50));
    c=colorbar;
    ylabel(c, '\rm \lambda','FontSize',16)
    ylabel('\it v','FontSize',16); title('Evolution of Simpson''s diversity index','FontSize',16);
    xlim([(years-1)*tf+1 years*tf]); 
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    subplot('Position',botpanelpos)
    i2 = imagesc(t((years-1)*tf+1:years*tf),1-y(:,1),G.Simp((years-1)*tf+1:years*tf)');set(gca, 'YDir', 'normal'); c=colorbar;
    ylabel(c,'\rm \lambda','FontSize',16);
    ylabel('\it \phi_L','FontSize',16); 
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')

    figure(8)
    clf;
    subplot('Position',toppanelpos)
    i1 = imagesc(t((years-1)*tf+1:years*tf),D.v(:,1),D.Shan((years-1)*tf+1:years*tf)'); set(gca, 'YDir', 'normal');
    colormap(jet(50));
    c=colorbar;
    ylabel(c, '\it H','FontSize',16)
    ylabel('\it v','FontSize',16); title('Evolution of Shannon''s diversity index','FontSize',16)
    xlim([(years-1)*tf+1 years*tf]); 
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    subplot('Position',botpanelpos)
    i2 = imagesc(t((years-1)*tf+1:years*tf),1-y(:,1),G.Shan((years-1)*tf+1:years*tf)');set(gca, 'YDir', 'normal'); c=colorbar;
    ylabel(c,'\it H','FontSize',16);
    ylabel('\it \phi_L','FontSize',16); 
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')
end 
