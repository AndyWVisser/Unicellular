function [simpD,simpG,shanD,shanG,Totsimp,Totshan]=D_G_diversity(D,G,t,tf,years)
    %This function calculates Simpson and Shannon index for one community of
    %diatoms, one community and fore one community composed of both diatoms
    %and generalists
    
    %Input : Diatoms and generalists structure, and time parameters 

    %preallocation
    D.BB = D.Bio; G.BB = D.BB; D.Psi = D.BB; G.Psi = D.BB; TotBB=D.BB; 
    simpD=zeros(tf*years,1); simpG = simpD; shanD=simpD; shanG = simpG; 
    Totsimp=zeros(tf*years,1); Totshan=Totsimp;

    for it=1:length(t)
        % diversity index calculation
        hd = squeeze(D.Bio(it,:,:));hg = squeeze(G.Bio(it,:,:));
        D.BB(it,:,:) = sum(hd(:)); G.BB(it,:,:) = sum(hg(:)); %total of all biomass at t=it
        TotBB(it,:,:) = D.BB(it,:,:)+G.BB(it,:,:);
        D.Psi(it,:,:) = D.Bio(it,:,:)./D.BB(it,:,:); G.Psi(it,:,:) = G.Bio(it,:,:)./G.BB(it,:,:); %proprotion of type i cells at t=it 

        simpD(it) = 1/sum(sum(D.Psi(it,:,:).^2)); simpG(it)= 1/sum(sum(G.Psi(it,:,:).^2)); 
        shanD(it) = exp(- sum(sum(D.Psi(it,:,:).*log(D.Psi(it,:,:))))); shanG(it) = exp(- sum(sum(G.Psi(it,:,:).*log(G.Psi(it,:,:)))));
        %Totsimp(it) = 1/sum(sum((D.Bio(it,:,:)./TotBB(it,:,:)).^2 + (G.Bio(it,:,:)./TotBB(it,:,:)).^2));
        Totsimp(it) = 1/sum(sum((D.Bio(it,:,:)./TotBB(it,:,:)).^2 + (G.Bio(it,:,:)./TotBB(it,:,:)).^2 ));
        %Totshan(it) =exp( - sum(sum((D.Bio(it,:,:)./TotBB(it,:,:)).*log(D.Bio(it,:,:)./TotBB(it,:,:)) + (G.Bio(it,:,:)./TotBB(it,:,:)).*log(G.Bio(it,:,:)./TotBB(it,:,:)))));
        Totshan(it) = exp(- sum(sum(D.Bio(it,:,:)./TotBB(it,:,:).*log(D.Bio(it,:,:)./TotBB(it,:,:))+G.Bio(it,:,:)./TotBB(it,:,:).*log(G.Bio(it,:,:)./TotBB(it,:,:)))));
    end 
    
end 