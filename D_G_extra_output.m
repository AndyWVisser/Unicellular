function [TTE_2, TTE_3, TTE_np,divmixo,divphoto]= D_G_extra_output(t,D,G,T,N,Si,tf,years) 

%This function starts like the ODE one, but only calculate variables useful
% for the calculation of TTE or productivity  (and diversity of phototrophs
% or mixotrophs)


global Imax Jmax x xo y 
global alphaF palat mp0 m2 Z 
global gmax P0 reminN reminS D0
global JR betaL betaN betaSi betaF
global L N0 S0 dil 

%preallocation
growthmixo=zeros(years*tf,Jmax,Imax); growthphoto=growthmixo; growthD=growthphoto; growthG=growthphoto;
G.mixo=zeros(years*tf,Jmax,Imax); G.photo = G.mixo;
mp=zeros(tf*years,Jmax,Imax); 
TTE_2=zeros(years*tf,1); TTE_3 = TTE_2; TTE_np = TTE_3;
divmixo = zeros(years*tf,1); divphoto = divmixo ; 
mixoTot = zeros(years*tf,Jmax,Imax); phoTot= mixoTot; 
for it=1:tf*years
    %===Reshaping biomass into a matrix for calculation===
     D.b= squeeze(D.Bio(it,:,:));
     G.b = squeeze(G.Bio(it,:,:));
     Nday=N(it);
     Siday=Si(it);
     Lday=L(it);
     Zday=Z(it);
    %===Trophic interactions===
     Fd = zeros(1,Imax); Fg = zeros(1,Imax);
     G.JF= zeros(Jmax,Imax);
     D.m = zeros(Jmax,Imax);     
     Bpred = sum(G.b.*G.phif);   %the biomass of predators [mg C / m-3] tin size classes r(i)
     for i = 1:Imax  % calculated once
        Fd(i) = sum(T(:,:,i).*D.b,'all');%food available from diatoms [mg C / m-3] to all generalists in size class x(i)
        Fg(i) = sum(T(1,:,i).*G.b,'all'); %food available from generalists [mg C / m-3] to all generalists in size class x(i)
        G.JF(:,i) = alphaF*G.phif(:,i)*(palat*Fd(i) + Fg(i));% J_F=A_F*F in size class x(i)

     end

      %===Uptake rates===
      G.JL = Lday.*G.AL; G.JN=Nday*G.AN;

      D.JL = Lday.*D.AL;  D.JN=Nday*D.AN; D.JSi= Siday*D.AN;

     %===Downregulation===
     D.Jresp = JR.*x; %Basal metabolism
     G.Jresp = JR.*x;

     D.limind = ones(Jmax,Imax);
     D.epsL= ones(Jmax,Imax);D.epsN = zeros(Jmax,Imax);D.epsSi = D.epsN;
     G.eps = max(0,min(1,(G.JL*(1-betaL)-G.Jresp-betaF*G.JF)./((betaN+G.CN).*G.JN)));

     %Different eventual cases :
     D.epsN1 = max(0,min(((1-betaL)*D.JL-D.Jresp)./((betaN + D.CN.*(1+betaSi./D.CS)).*D.JN),1)); %C limiting
     D.epsSi1 = max(0,min(((1-betaL)*D.JL-D.Jresp)./((betaSi + D.CS + betaN.*D.CS./D.CN).*D.JSi),1));

     D.epsSi2 = max(0,min(D.CN.*D.JN./(D.JSi.*D.CS),1)); %N limiting 
     D.epsL2 = max(0,(D.Jresp + D.CN.*D.JN +betaN.*D.JN + betaSi*D.CN.*D.JN./D.CS)./((1-betaL)*D.JL));

     D.epsN3 = max(0,min(1,D.CS.*D.JSi./(D.CN.*D.JN)));%Si limiting
     D.epsL3 =max(0,(D.Jresp + (D.CS+betaSi).*D.JSi + betaN*D.CS.*D.JSi./D.CN)./((1-betaL)*D.JL));    

     for i=1:Jmax*Imax
        if D.epsN1(i)<1 && D.epsSi1(i)<1 %carbon is limiting
            D.epsN(i) = D.epsN1(i);
            D.epsSi(i) = D.epsSi1(i);
        elseif D.epsSi2(i)<1 && D.epsL2(i) <1 %Nitrogen is limiting
            D.epsSi(i) = D.epsSi2(i); 
            D.limind(i) = 2 ; 
            D.epsN(i) = 1; 
        elseif D.epsN3(i)<1 && D.epsL3(i)<1 %Silica is limiting 
            D.epsN(i) = D.epsN3(i);
            D.epsSi(i) = 1; 
            D.limind(i) = 3;      
        end           
     end 


    D.JNreg = D.epsN.*D.JN;
    D.JSireg = D.epsSi.*D.JSi;
    D.Ccons = max(0,(1-betaL)*D.JL-D.Jresp-betaL*D.epsL.*D.JL-betaN*D.JNreg-betaSi*D.JSireg); 
    D.Ncons = D.CN.*D.JNreg;
    D.Sicons = D.CS.*D.JSireg;
    D.Jeff = min(min(D.Ncons,D.Sicons),D.Ccons); 
    D.Jeff = max(D.Jeff,0); %Effective uptake (cannot be negative)


    G.JNreg = G.JN.*G.eps ; 
    G.Ccons = max(0,G.JF+(1-betaL)*G.JL-G.Jresp-betaF*G.JF-betaN*G.JNreg);
    G.Ncons = G.CN.*G.JNreg+G.JF ;
    G.Jeff = min(G.Ccons,G.Ncons);
    
    %==== Cells division rate ====%
    D.g = gmax.*D.Jeff./(D.Jeff+gmax.*ones(Jmax,Imax));
    G.g = gmax.*G.Jeff./(G.Jeff+gmax.*ones(Jmax,Imax));

    %==== Grazing by copepods ====%
     D.mp = zeros(Jmax,Imax); G.mp = zeros(Jmax,Imax);
        for i=1:Jmax*Imax
            if x(i)> 10^(-2)        G.mixo=zeros(tf,Jmax,Imax); G.photo = G.mixo;
        mp=zeros(tf,Jmax,Imax); 
                D.mp(i)= mp0*(1-palat)*Zday*D.r(i).^(-3/4);
                G.mp(i)= mp0*Zday*G.r(i).^(-3/4);
            end 
        end 
    
    for i = 1:Imax  
        D.m = D.m + alphaF*Bpred(i).*T(:,:,i).*(ones(Jmax,Imax)-D.g(:,i)); %predation by larger cells for diatoms 
    end

      G.m = D.m(1,:);
      D.m = palat*D.m; % mortality rates [day-1] due to internal predation 
      D.m(D.m<0)=0;G.m(G.m<0)=0;

    
    %Growthrates for mixotrophs, phototrophs, diatomos & generalists
    growthmixo(it,:,:)= G.g.*G.phif; %mixotrophs
    growthphoto(it,:,:)= D.g + G.g.*(ones(Jmax,Imax)-G.phif);%phototrophs
    growthD(it,:,:)=D.g;
    growthG(it,:,:)=G.g;
    
    %Reshaping 
     G.mixo(it,:,:)=reshape(squeeze(G.Bio(it,:,:)).*(ones(Jmax,Imax)-G.phil),[1 Jmax Imax]); %biomass of mixotroph generalists
     G.photo(it,:,:)= reshape(squeeze(G.Bio(it,:,:)).*G.phif,[1 Jmax Imax]); %biomass of phototroph generalists
    mp(it,:,:)= reshape(D.mp + G.mp,[1 Jmax Imax]); %copepods predation rate
    
    
    %Trophic Transfer Efficiency 
    TTE_2(it)=100*sum(squeeze(G.mixo(it,:,:)).*G.JL,'all')./sum(squeeze(G.photo(it,:,:)+D.Bio(it,:,:)).*D.JL,'all'); %1st level: phototrophs to mixotrophs
    %TTE_2(it)=100*sum(growthmixo(it,:,:),'all')/sum(growthphoto(it,:,:),'all');
    %TTE_2(it)=100*sum(squeeze(G.mixo(it,:,:)).*G.g./(squeeze(G.photo(it,:,:)+D.Bio(it,:,:)).*D.g),'all');
    TTE_3(it)= 100*sum(mp(it,:,:),'all')/sum(growthmixo(it,:,:),'all'); %2nd level : mixotrophs to copepods
    TTE_np(it) = 100*sum(mp(it,:,:),'all')/sum(growthphoto(it,:,:),'all'); %non progressive, only one level : phototrophs to copepods
    
    
    
    %Diversity of mixotrophs and of phototrophs
    hmix= squeeze(G.mixo(it,:,:)); mixoTot(it,:,:) = sum(hmix(:));
    hpho = squeeze(G.photo(it,:,:)+D.Bio(it,:,:)); phoTot(it,:,:)= sum(hpho(:));
    
    divmixo(it) = exp( - sum(sum(G.mixo(it,:,:)./mixoTot(it,:,:).*log(G.mixo(it,:,:)./mixoTot(it,:,:))))); 
    divphoto (it)=exp( - sum(sum((G.photo(it,:,:)+D.Bio(it,:,:))./phoTot(it,:,:).*log((G.photo(it,:,:)+D.Bio(it,:,:))./phoTot(it,:,:)))));
  
    
end


end 