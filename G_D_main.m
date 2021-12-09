%This code simulates runs for STEADY or SEASONAL environment 
%Things need to be defined beforehand : 
% - Number of size classes Imax, number of trophic strategy classes Jmax
% - steady (seasonal=0) or seasonal (seasonal=1) environmental conditions
% - run time in terms or years 

%Structure of the code
% 1. define state variables
% 2. define environmental conditions
% 3. solve the equations with ODE
% 4. use that output for calculation of interest variables (diversity,
% trophic transfer, etc)
% 5. plot the figures needed

 %==== Dialog box ====%
rep = inputdlg({'Size classes (number)','Trait classes (number)','seasonal (0/1)','years (number)'},'Requirements',[1 12;1 12;1 12;1 12]);

%The following global variables are needed for the ODE function
global Imax Jmax x xo y 
global alphaF palat mp0 m2 Z 
global gmax P0 reminN reminS D0
global JR betaL betaN betaSi betaF
global L N0 S0 dil 


%==== Parameters ====%
%Carbon,Nitrogen,Si density
cb = 0.11e-6; % [µgC/µm^3] carbon content of cytoplasm (from Sthratmann,1967)
cm = 0.6e-6;  % [µgC/µm^3] carbon content of membranes
gb = 7;       % C:N ratio of cytoplasm [µgC/µgN]
gm = 38;      % C:N ratio of membrane [µgC/µgN]
Sishell = 10.3E-6; % silicate concentration in diatom shell [µg Si µm-3]
%cL = 0.0318*(gb/cb)^(-2/3);   % light affinity coefficient [µgC/day/(W m^-2)/µm^2]

%Thickness
h = 0.008;    % cell membranes [µm]
tshell = 0.01;     % shell [µm] 

%Affinity coefficients
kappa = 0.2;     % light affinity investment return [µm^-1] (also called absorption coefficient)
%cN = 2.89E-4*(gb/cb)^(-1/3); % nutrient affinity coefficient [L/day/µm]
%cP = 0.227*(gb)^(-1); % food affinity coefficient [L/day/µgC]
%eP = 0.3408*(gb)^(-1) ; % nutrient affinity investment return [L/day/µgC]

alphaN = 0.68;  % L day-1 (µg C)-1 µm2  
alphaL = 0.18 * 1E-12 * 24 * 3600*1.875;   % quantum yield (scaled)[µgC day-1 (µE m-2 s-1)-1 µm-2] 
qyield = 0.07; % quantum yield {µg C / µmol photons PAR]
%useful for the calculation of nutrients affinity
rstar = 2;  % crossover from r3 tp r absolute nutrient affinity 
xstar = 4*pi/3*rstar^3*cb; % mass equivalence for rstar

%Food 
alphaF = 0.018; %0.018 food affinity coefficient[L d-1 (µgC)-1]
palat = 0.6; % diatom palatibility factor
epsilon = 7; % predator-prey size ratio (linear dimension)
width = 1; % in orders of magnitude
mp0 = 0.1; % grazing constant 
m2 = 0.01; % quadratic viral lysis [(µC L^-1)^-1 d^-1]
gmax = 1; %maximal growth rate [d^-1]
Z0 = 0.06; %copepods [µgC L^-1] (predation by higher trophic levels)
reminS = 0.8; reminN = 0.1; %remineralization rates 

%Metabolic costs 
JR = 0.04; % basal respiration [µgC µgC^-1 day^-1] (beta0)
betaL = 0.08;   % cost of light harvesting  [µgC µgC^-1]
betaN = 0.45;   % cost of nitrogen uptake  [µgC µgN^-1]
betaSi = 0.45;  % cost of silicate uptake  [µgC µgSi^-1]
betaF = 0.45;   % cost of particulate matter uptake  [µgC µgC^-1]

%Classes
Imax = str2num(rep{1}); %Number of size classes 
Jmax = str2num(rep{2}); %Number of classes for vacuolation or trophic strategies 


%==== Primary trait space ====%
xo = 10.^linspace(-7,1,Imax); %cell carbon mass ranging from 10^-7 to 10 µgC cell^-1
[x,y] = meshgrid(xo,linspace(0,1,Jmax)); %y is phagotrophic feeding for generalists
vo = linspace(0,1,Jmax+1); % vacuole size for diatoms 
[x,v] = meshgrid(xo,vo(1:Jmax));
r = (3/(4*pi*cb)*x./(1-v)).^(1/3); %celle equivalent spherical radius [µm]

phib = 0; %Investment in biomass photosynthesis
xm = 4*pi*r.^2*h*cm.*(1 + v.^(2/3)); %Carbon mass in membranes
phim = xm./x; %Carbon amount apportioned to membrane mass

%Generalists
G.x = x; G.r = ones(Jmax,1)*r(1,:); %(only radius with v=0)
G.phim = ones(Jmax,1)*phim(1,:); %phi_M(x,0) (first line of phim corresponds to v=0) 
phil = 1 - phim - phib - y; phil(phil<0) = 0; 
phif = y; phif(phil<0) = 0;
G.phil = phil; G.phif = phif;

%Diatoms
D.x=x;D.v=v;D.r=r; 
D.phim=phim;
D.phim(D.phim > 1) = 1;
phil = 1 - phim - phib; phil(phil<0) = 0;
D.phil = phil;

%N and Si requirements
D.shell = 2.28E-3*(D.r).^(0.72); % thickness of diatom shell (t=1.62V^0.24)
D.s = 4*pi*Sishell*(D.r.^2).*D.shell; % diatom Si requirement [µg Si]
D.n = xm/gm + (D.x - D.phim.*D.x)/gb; % diatom N requirement [µg N]
%D.n(D.n<0)=0;
G.n = G.phim.*G.x/gm + (G.x-G.phim.*G.x)/gb; %generalist N requirement [µg N]

%C:Si ratio
D.CS = D.x./D.s; 
%C:N ratio 
D.CN= D.x./D.n;
G.CN = G.x./G.n; 

%==== Affinitites ====%
D.AL = alphaL*pi*(D.r.^2).*(1 - exp(-kappa*D.phil.*D.r.*(1-D.v)))./D.x; D.AL = D.AL*0.9; %0.9 factor is an attempt of correction to match Mathilde's value
D.AN = 4/3*pi*cb*alphaN*(D.r)./(1 + (D.x/xstar).^(-2/3))./D.x; D.AN = D.AN;
D.AN(D.phim == 1) = 0;

G.AL = alphaL*pi*(G.r.^2).*(1 - exp(-kappa*G.phil.*G.r))./G.x;  
G.AN = alphaN*((G.r).^(-2))./(1 + (G.r/rstar).^(-2)); 
G.AF = alphaF*G.phif; 

G.JF = zeros(size(G.AF));

%==== Trophic interactions ====%
Theta = @(rpredator,rprey) exp(-(log(epsilon*rprey./rpredator)).^2/(2*width^2)); %represents size preference for smaller cells 

T = zeros(Jmax,Imax,Imax);  
for i = 1:Imax  % calculated once 
    T(:,:,i) = Theta(r(1,i),D.r); %3D matrix with all possible predator prey interactions 
end

%==== Time parameters ====%
tf = 365;
tspan=[1:tf];
years = str2num(rep{4});
tspan_y = [1:tf*years];

%==== Steady or seasonal environmental forcing ====%
seasonal = str2num(rep{3});

%==== Initial conditions ====%
%Plankton
Binit = 1;  % Initial total biomass abundance [mgC m^-3]
BGtoD = 0.7;  % biomass split BGtoD in generalists, 1- BGtoD in diatoms
P0 = ones([Jmax Imax])*0.09;% %Background plankton biomass concentration [mgC L^-1] 
D0 = P0; %backroung concentration for diatoms 

%little adjusment
for i=1:Imax*Jmax
    if D.AN(i)==0
        D0(i)=0; %taking out the "small cell size, big vacuole" diatoms that should'nt exist
    end 
end

D.b = ones(size(phil))*Binit*(1-BGtoD)./(Imax*Jmax); %initial total biomass abundance for diatoms 
G.b = ones(size(phil))*Binit* BGtoD./(Imax*Jmax); %initial total biomass abundance for generalists

%correct format for the ODE solver
B0 = [reshape(D.b,[Jmax*Imax 1]);reshape(G.b,[Jmax*Imax 1])]; %initial vector for the ODE
maxind = 2*Imax*Jmax; %size of the initial vector

%Environmental conditions
lat = 55; %possible to change latitude
[L,dil,Z,N0,S0]=D_G_physical_settings(lat,years,seasonal,tspan,tf); 

%==== Solving the equations ====%
if seasonal ==0
    D.allB = []; G.allB = D.allB; 
    D.Bio = zeros(tf*years,Jmax,Imax); G.Bio = D.Bio;
    Ls = L; dils= dil; Zs = Z; N0s = N0; S0s= S0; 
    for s=1:3 %for the 3 steady states
        L=Ls(s); dil = dils(s); Z = Zs(s); N0=N0s(s); S0 = S0s(s);
        
         options = odeset('NonNegative',[1:maxind],'RelTol',1e-4,'AbsTol',1e-4); 
        [t, B] = ode45(@(t,B) D_G_ode_func(t,B,D,G,T,seasonal),tspan_y,B0,options);
        [No_use,HD,HG] = D_G_ode_func(t,B,D,G,T,seasonal)  ; %output of D and G
        for it=1:length(t)
            %reshaping biomass into a 3D matrix
            D.Bio(it,:,:) = reshape(squeeze(B(it,1:Jmax*Imax)),[Jmax Imax]); 
            G.Bio(it,:,:) = reshape(squeeze(B(it,Jmax*Imax+1:Jmax*2*Imax)),[Jmax Imax]);
        end 
        D.allB = [D.allB; D.Bio]; G.allB =[G.allB; G.Bio]; %Biomass for the 3 steady states 
    end 

else 
    B0 = [B0;N0;S0]; %if seasonal, nutrients concentration also varies over time -> (2*Jmax*Imax + 2) vector : nutrients variables at the end 
    maxind = Jmax*2*Imax+2; %size of the vector
    
    options = odeset('NonNegative',[1:maxind],'RelTol',1e-4,'AbsTol',1e-4); 
    [t, B] = ode45(@(t,B) D_G_ode_func(t,B,D,G,T,seasonal),tspan_y,B0,options);
    [No_use,HD,HG] = D_G_ode_func(t,B,D,G,T,seasonal)  ; %output of D and G

    D.Bio = zeros(tf*years,Jmax,Imax); G.Bio = D.Bio; G.photo=G.Bio; G.mixo=G.photo;
    for it=1:length(t)
        %reshaping biomass into a 3D matrix
        D.Bio(it,:,:) = reshape(squeeze(B(it,1:Jmax*Imax)),[Jmax Imax]); 
        G.Bio(it,:,:) = reshape(squeeze(B(it,Jmax*Imax+1:Jmax*2*Imax)),[Jmax Imax]);
        G.photo(it,:,:) = squeeze(G.Bio(it,:,:)).*G.phil;
        G.mixo(it,:,:) = squeeze(G.Bio(it,:,:)).*(ones(Jmax,Imax)-G.phil);
    end

    %nutrients concentration
    N = B(:,2*Jmax*Imax+1);
    Si = B(:,2*Jmax*Imax+2);
    
    %diversity indexes
    [D.Simpson, G.Simpson, D.Shannon,G.Shannon, Totsimpson, Totshannon]=D_G_diversity(D,G,t,tf,years);

    
end   


    D.Btot = squeeze(sum(D.Bio,2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.Bio,2)); %sum for all trophic strategies for one size class
    D.Btotot = sum(D.Btot,2); %total biomass for diatoms at each time step (size = [years*365 1])
    G.Btotot = sum(G.Btot,2);
    G.photo_tot = sum(squeeze(sum(G.photo,2)),2);
    G.mixo_tot = sum(squeeze(sum(G.mixo,2)),2);    
    
%==== Trophic transfer ====%
    if seasonal==1
        [TTE_2, TTE_3, TTE_np,divmixo,divphoto]= D_G_extra_output(t,D,G,T,N,Si,tf,years);
        
                %==== Productivity ====%
        D.prod = zeros(years*tf,Jmax,Imax);
        G.prod = zeros(years*tf,Jmax,Imax);
        for it=1:years*tf
            D.prod(it,:,:)=L(it)*D.AL.*reshape(D.Bio(it,:,:),[Jmax Imax]);
            G.prod(it,:,:)=L(it)*G.AL.*reshape(G.Bio(it,:,:),[Jmax Imax]);
        end

        D.prod = sum(squeeze(sum(D.prod,2)),2)*10^(-3);
        G.prod = sum(squeeze(sum(G.prod,2)),2)*10^(-3);

    end 

%==== Plot ====%

D_G_plot; 
