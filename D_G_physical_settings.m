function [L,dil,Z,N0,S0]=D_G_physical_settings(lat,years,seasonal,tspan,tf)
    %This function computes the environmental conditions, according to the
    %type of simulation (steady/seasonal), the latitude and time parameters

if seasonal==1
    %initial nutrient concentration
    N0 = 1200;  % background N concentration ([µgSi L-1] ) (also initial N concentration in the photic layer) 
    S0 = 1000; % background Si concentration  (also initial Si concentration in the photic layer) 
  
    
    % light L [µmol photons m^-2 s^-1]
    clouds=1-(6/8); %attenuation by clouds
    I=daily_insolation(0,lat,1:365,1); %surface light
    L=I.*clouds.*0.4.*4.6; %add the cloud cover and convert to the right unit
    L=L';
    L=[repmat(L(1:end),years,1)];
    
    %dilution rate [d-1]
    dil1 = 105*(1-cos(2*pi*tspan/365))+25; dil1=dil1'; 
    dil = 0.9*(1+(-dil1+min(dil1))/210)+0.1;
    dil = repmat(dil(1:end),years,1);
    
    %grazers abundance [µgC L-1] 
    Z=ones(tf,1)*5; Z(1:75)=5; Z(76:135)=0.25*tspan(76:135)-13.75; Z(136:239)=5/104*tspan(136:239)+13.5; 
    Z(240:269)=-0.6*tspan(240:269)+168.4; Z(270:365)=-2/51*tspan(270:365)+5+2*320/51;       
    Z=smooth(Z,50);
    Z= repmat(Z(1:end),years,1);

else % 3 steady states (winter, spring, summer)
     L = [50 120 220]; 
    dil= [0.5 0.2 0.1]; 
    Z = [5 10 20]; 
    N0 = [50 30 8];  % background N concentration  [µgN L-1] (not sure about units) 
    S0 = [80 20 8];  % background Si concentration  [µgSi L-1] 
    
end 

end 