%This script plot all the necessary figures, for steady or seasonal simulation:
%Dialog box asks if you're interest by diversity plot and if so which diversity index
%would you like
toppanelpos = [.13,.51,.78,.4];
botpanelpos = [.13,.09,.78,.4];

%for axes display
pow = char('10^{-7}','10^{-5}','10^{-3}','10^{-1}','10');
power = char('10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10');

%==== Steady state ====% 
if seasonal==0 
    pos1= [.06,.66,.27,.27]; pos2= [.36,.66,.27,.27]; pos3= [.66,.66,.27,.27];
    pos4= [.06,.36,.27,.27]; pos5= [.36,.36,.27,.27]; pos6= [.66,.36,.27,.27];
    pos7= [.06,.06,.27,.27]; pos8= [.36,.06,.27,.27]; pos9= [.66,.06,.27,.27];
  
    figure(20)
    clf;
    set(gcf,'color','w');
    subplot('Position',pos1)
    D.Btot = squeeze(sum(D.allB(1:365,:,:),2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.allB(1:365,:,:),2));
    D.Btotot = sum(D.Btot(end,:),2); G.Btotot = sum(G.Btot(end,:),2);
    plot(log10(xo),D.Btot(end,:),log10(xo),G.Btot(end,:)); 
    legend('Diatoms','Generalists');
    title('Winter'); ylabel('Biomass [µgC L^{-1}]');
    xlim([-7 1]);xticks([-7 -5 -3 -1 1]); xticklabels(pow); 

    
    subplot('Position',pos4)
    futG=reshape(G.allB(365,:,:),[Jmax Imax]);
    imagesc(log10(xo),1-y(:,1),futG);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    ylabel(c, 'P \rm [µgC L^{-1} ]')
    ylabel('\rm \phi_L'); set(gca,'XTick', []);
   
    
    subplot('Position',pos7)
    futD=reshape(D.allB(365,:,:),[Jmax Imax]);
    imagesc(log10(xo),D.v(:,1),futD);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    ylabel('\rm v'); xlabel('\it x \rm [µg C]'); xticks([-7 -5 -3 -1 1]); xticklabels(pow);
    
    subplot('Position',pos2)
    D.Btot = squeeze(sum(D.allB(366:730,:,:),2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.allB(366:730,:,:),2));
    D.Btotot = [D.Btotot sum(D.Btot(end,:),2)]; G.Btotot =[G.Btotot sum(G.Btot(end,:),2)]; 
    plot(log10(xo),D.Btot(end,:),log10(xo),G.Btot(end,:)); title('Spring');xlim([-7 1]); xticks([-7 -5 -3 -1 1]); xticklabels(pow);
    
    subplot('Position',pos5)
    futG=reshape(G.allB(730,:,:),[Jmax Imax]);
    imagesc(log10(xo),1-y(:,1),futG);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    axis off
    
    subplot('Position',pos8)
    futD=reshape(D.allB(730,:,:),[Jmax Imax]);
    imagesc(log10(xo),D.v(:,1),futD);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    xlabel('\it x \rm [µg C]'); xticks([-7 -5 -3 -1 1]); xticklabels(pow); set(gca,'YTick', []);
   
    
    subplot('Position',pos3)
    D.Btot = squeeze(sum(D.allB(731:1095,:,:),2)); %sum for all vacuoles sizes for one size class
    G.Btot = squeeze(sum(G.allB(731:1095,:,:),2));
    D.Btotot = [D.Btotot sum(D.Btot(end,:),2)]; G.Btotot =[G.Btotot sum(G.Btot(end,:),2)]; 
    plot(log10(xo),D.Btot(end,:),log10(xo),G.Btot(end,:)); title('Summer');xlim([-7 1]); xticks([-7 -5 -3 -1 1]); xticklabels(pow);
    
    subplot('Position',pos6)
    futG=reshape(G.allB(1095,:,:),[Jmax Imax]);
    imagesc(log10(xo),1-y(:,1),futG);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    axis off
    
    subplot('Position',pos9)
    futD=reshape(D.allB(1095,:,:),[Jmax Imax]);
    imagesc(log10(xo),D.v(:,1),futD);
    set(gca, 'YDir', 'normal');colormap(jet(50)); c=colorbar;
    xlabel('\it x \rm [µg C]'); xticks([-7 -5 -3 -1 1]); xticklabels(pow); set(gca,'YTick', []);
   
% 
%     figure(2)
%     clf;
%     set(gcf,'color','w');
%     plot(D.allSimp,D.Btotot,G.allSimp,G.Btotot); 
%     xlabel('D'); ylabel('Total biomass [mgC m^{-3}]');
%     legend('diatoms','generalists')
%   
    
    
%==== Seasonal cycle ====%    
else 
    pos1=[.06,.06,.28,.8] ; pos2=[.38,.06,.28,.8]  ; pos3= [.7,.06,.28,.8] ; 
    
    %Dialog box 
    prompt ={'Do you want the diversity index ? (yes/no)','Which one ? (Simpson/Shannon'} ; 
    dlgtitle= 'Diversity';
    dims = [1 40;1 40];
    answer = inputdlg(prompt,dlgtitle,dims);
    user_v = string(answer{1});
    user_w = string(answer{2});
    
    D.Btot_y=D.Btot((years-1)*tf+1:years*tf,:);
    G.Btot_y=G.Btot((years-1)*tf+1:years*tf,:);
    Bmax=max(max(max(G.Btot_y)),max(max(D.Btot_y)));
    
    %Biomass evolution for diatoms
    figure(1)
    set(gcf,'color','w');
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),D.Btot_y',[0 Bmax]); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf])
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    A=jet(120);
    A(1,:)=1;
    colormap(A);
    C=colorbar;
    xlabel('Months')
    xlabel(C,'Biomass \rm [µgC L^{-1}]');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp   
    yticks([-7 -6 -5 -4 -3 -2 -1 0 1]); yticklabels(power); 

     %Biomass evolution for generalists
    figure(2)
    set(gcf,'color','w');
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),G.Btot_y',[0 Bmax]); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf])
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    colormap(A);
    C=colorbar;
    xlabel('Months')
    xlabel(C,'Biomass \rm [µgC L^{-1}]');
    ylabel('Mass \rm [µgC]');
    set(gca,'FontSize',14);
    shading interp   
    yticks([-7 -6 -5 -4 -3 -2 -1 0 1]); yticklabels(power); 
    
    %Evolution of nutrients concentration
    figure(3)
    set(gcf,'color','w');
    plot(t((years-1)*tf+1:years*tf),N((years-1)*tf+1:years*tf),'b',t((years-1)*tf+1:years*tf),Si((years-1)*tf+1:years*tf),'g','LineWidth',2)
    legend('N \rm [µgN L^{-1}]','Si \rm [µgS L^{-1}]');
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')
    

    
    
%plotting dominant traits 
    D.Btrait = zeros(tf*years,Imax);
    for it=1:tf*years
        for i=1:Imax
            k=max(D.Bio(it,:,i));
            kbis= max(find(D.Bio(it,:,i)==k));
            D.Btrait(it,i) = D.v(kbis,1);
            if D.Btot(it,i)< (Bmax/120) %below, biomass considered non-existent
                D.Btrait(it,i) = NaN;
            end
        end 
    end 

    G.Btrait = zeros(tf*years,Imax);
    for it=1:tf*years
        for i=1:Imax
            k=max(G.Bio(it,:,i));
            kbis= max(find(G.Bio(it,:,i)==k));
            G.Btrait(it,i) = 1-y(kbis,1); 
             if G.Btot(it,i)< (Bmax/120)
                G.Btrait(it,i) = NaN;
            end
        end 
    end 
    
    %Dominant trait values
    figure(4)
    clf;
    set(gcf,'color','w');
    D.Btrait_y=D.Btrait((years-1)*tf+1:years*tf,:);
    imAlpha=ones(Imax,tf);
    imAlpha(isnan(D.Btrait_y'))=0;
    set(gca,'color',0*[1 1 1]);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),D.Btrait_y','AlphaData',imAlpha); 
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf]);
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    futcmap = parula(120);
    colormap(futcmap)
    C=colorbar;
    xlabel('Months')
    xlabel(C,'\it v (%)');
    ylabel('Mass \rm [µgC]'); yticks([-7 -6 -5 -4 -3 -2 -1 0 1]); yticklabels(power);
    set(gca,'FontSize',14);
    shading interp 

    figure(5)
    clf;
    set(gcf,'color','w');
    G.Btrait_y=G.Btrait((years-1)*tf+1:years*tf,:);
    imAlpha=ones(Imax,tf);
    imAlpha(isnan(G.Btrait_y'))=0;
    set(gca,'color',0*[1 1 1]);
    ie = imagesc(t((years-1)*tf+1:years*tf),log10(xo),G.Btrait_y','AlphaData',imAlpha);
    set(gca, 'YDir', 'normal');
    xlim([(years-1)*tf+1 years*tf]);
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr); %x axis : time in months
    xticklabels(monthstr);
    C=colorbar;
    xlabel('Months')
    xlabel(C,'\it \phi_L');
    ylabel('Mass \rm [µgC]'); yticks([-7 -6 -5 -4 -3 -2 -1 0 1]); yticklabels(power);
    set(gca,'FontSize',14);
    shading interp 
    
    %plotting the evolution of the total biomass over a seasonal cycle 
    figure(6)
    set(gcf,'color','w');
    plot(t((years-1)*tf+1:years*tf),D.Btotot((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),G.Btotot((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),G.mixo_tot((years-1)*tf+1:years*tf),'--',t((years-1)*tf+1:years*tf),G.photo_tot((years-1)*tf+1:years*tf),'-.','LineWidth',2); 
    legend('diatoms','generalists','mixotroph generalists','phototroph generalists');
    ylabel('Biomass [µgC L^{-1}]');
    xlim([(years-1)*tf+1 years*tf]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')
    
    %plotting trophic transfer efficiency over a seasonal cycle
    figure(15)
    set(gcf,'color','w');
    plot(t(1:tf),TTE_2((years-1)*tf+1:years*tf),t(1:365),TTE_3((years-1)*tf+1:years*tf),t(1:365),TTE_np((years-1)*tf+1:years*tf),'LineWidth',2);
    legend('TTE2','TTE3','TTE not progressive');
    ylabel('[%]');
    xlim([1 365]);   
    monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335];
    xticks(monthnr) 
    xticklabels(monthstr)
    xlabel('Months')
    
 
%     figure(16)
%     clf;
%     set(gcf,'color','w');
%     plot(Totshannon((years-1)*tf+1:years*tf),TTE_2((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),TTE_3((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),TTE_np((years-1)*tf+1:years*tf),'LineWidth',2);
%     %title('Trophic transfer efficiency as a function of diversity');
%     legend('TTE2','TTE3','TTE not progressive'); 
%     xlabel('Diversity \it E(H)');
%     ylabel('[%]');
%     xlim([min(Totshannon) max(Totshannon)]);
    
    %Trophic transfer efficiency as a function of diversity
    figure(16)
    clf;
    set(gcf,'color','w');
    subplot('Position',pos1);
    plot(divphoto((years-1)*tf+1:years*tf),TTE_2((years-1)*tf+1:years*tf),divmixo((years-1)*tf+1:years*tf),TTE_2((years-1)*tf+1:years*tf),'LineWidth',2);
    legend('mixotrophs diversity','phototrophs diversity')
    ylabel('Trophic Transfer Efficiency [%]');
    xlabel('Diversity \it E(H)');
    
    subplot('Position',pos2);
    plot(divmixo((years-1)*tf+1:years*tf),TTE_3((years-1)*tf+1:years*tf),'LineWidth',2);
    xlabel('Diversity \it E(H)');

    subplot('Position',pos3)
    plot(divphoto((years-1)*tf+1:years*tf),TTE_np((years-1)*tf+1:years*tf),'LineWidth',2);
    xlabel('Diversity \it E(H)');
    
% Phototrophs and mixotrophs diversity
%     figure(7)
%             clf;
%             set(gcf,'color','w');
%             plot(t((years-1)*tf+1:years*tf),divphoto((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),divmixo((years-1)*tf+1:years*tf),'LineWidth',2);
%             legend('photo','mixo'); 
%             ylabel('\it E(H)','FontSize',16); %title('Evolution of Shannon''s diversity index','FontSize',16)
%             xlim([(years-1)*tf+1 years*tf]); 
%             monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
%             monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
%             xticks(monthnr) 
%             xticklabels(monthstr)
%             xlabel('Months');

    
    %figure(17)
    % clf;
% filename = 'prodG.gif';
% for tt=1:365
% 
%      plot(log10(xo),squeeze(G.prod(tt+365,:,:)))
%      drawnow
%      frame = getframe(2);
%      im = frame2im(frame);
%      [imind,cm] = rgb2ind(im,256);
%      if tt == 1;
%            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%      else
%            imwrite(imind,cm,filename,'gif','WriteMode','append');
%      end
%  end 

%Productivity as a function of diversity
figure(19)
clf;
 set(gcf,'color','w');
 plot(D.Shannon((years-1)*tf+1:years*tf),D.prod((years-1)*tf+1:years*tf),G.Shannon((years-1)*tf+1:years*tf),G.prod((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),D.prod((years-1)*tf+1:years*tf)+G.prod((years-1)*tf+1:years*tf),'LineWidth',2);
 legend('diatoms','generalists','all'); 
 xlabel('Diversity \it E(H)');
 ylabel('Q [mgC d^{-1} L^{-1}]');
 xlim([min(D.Shannon) max(Totshannon)]);
                
                
       
    
    if strcmp(user_v,'yes')
        if strcmp(user_w,'Simpson')
            %Plotting the diversity indexes             
            figure(7)
            clf;
            set(gcf,'color','w');
            plot(t((years-1)*tf+1:years*tf),D.Simpson((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),G.Simpson((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),Totsimpson((years-1)*tf+1:years*tf),'LineWidth',2);
            legend('diatoms','generalists','all'); 
            ylabel('\it D','FontSize',16); %title('Evolution of Simpson''s diversity index','FontSize',16)
            xlim([(years-1)*tf+1 years*tf]); 
            monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
            monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
            xticks(monthnr) 
            xticklabels(monthstr)
            xlabel('Months')
            
                %Simpson Biodiversity as function of physical forcings
                figure(8)
                clf;
                set(gcf,'color','w');
                subplot('Position',pos1);
                plot(L(1:365),D.Simpson((years-1)*tf+1:years*tf),'.',L((years-1)*tf+1:years*tf),G.Simpson((years-1)*tf+1:years*tf),'.',L((years-1)*tf+1:years*tf),Totsimpson((years-1)*tf+1:years*tf),'.');
                legend('diatoms','generalists','all'); 
                ylabel('Diversity \it D');
                xlabel('Light L [µmol photons m^{-2} s^{-1}]');
                xlim([min(L) max(L)]);

                subplot('Position',pos2);
                plot(dil((years-1)*tf+1:years*tf),D.Simpson((years-1)*tf+1:years*tf),'.',dil((years-1)*tf+1:years*tf),G.Simpson((years-1)*tf+1:years*tf),'.',dil((years-1)*tf+1:years*tf),Totsimpson((years-1)*tf+1:years*tf),'.');
                set(gca,'YTick', []);
                xlabel('Dilution rate [d^{-1}]');
                xlim([min(dil) max(dil)]);

                subplot('Position',pos3);
                plot(Z(1:365),D.Simpson((years-1)*tf+1:years*tf),'.',Z(1:365),G.Simpson((years-1)*tf+1:years*tf),'.',Z(1:365),Totsimpson((years-1)*tf+1:years*tf),'.');
                set(gca,'YTick', []);
                xlabel('Copepods abundance [µgC L^{-1}]');
                xlim([min(Z) max(Z)]);
                

            
        else %Shannon
            
            %Evolution of Shannon''s diversity index
            figure(7)
            clf;
            set(gcf,'color','w');
            plot(t((years-1)*tf+1:years*tf),D.Shannon((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),G.Shannon((years-1)*tf+1:years*tf),t((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),'LineWidth',2);
            legend('diatoms','generalists','all'); 
            ylabel('\it E(H)','FontSize',16); 
            xlim([(years-1)*tf+1 years*tf]); 
            monthstr = {'J','F','M','A','M','J','J','A','S','O','N','D'};
            monthnr  = [  1   32  60  91  121 152 182 213 244 274 305 335]+ones(1,12)*tf*(years-1);
            xticks(monthnr) 
            xticklabels(monthstr)
            xlabel('Months');
                
                %Shannon Biodiversity as function of physical forcings
                figure(20)
                clf;
                set(gcf,'color','w');
                subplot('Position',pos1);
                plot(L((years-1)*tf+1:years*tf),D.Shannon((years-1)*tf+1:years*tf),'.',L((years-1)*tf+1:years*tf),G.Shannon((years-1)*tf+1:years*tf),'.',L((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),'.');
                legend('diatoms','generalists','all'); 
                ylabel('Diversity \it E(H)');
                xlabel('Light L [µmol photons m^{-2} s^{-1}]');
                xlim([min(L) max(L)]);

                subplot('Position',pos2);
                plot(dil((years-1)*tf+1:years*tf),D.Shannon((years-1)*tf+1:years*tf),'.',dil((years-1)*tf+1:years*tf),G.Shannon((years-1)*tf+1:years*tf),'.',dil((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),'.');
                set(gca,'YTick', []);
                xlabel('Dilution rate [d^{-1}]');
                xlim([min(dil) max(dil)]);

                subplot('Position',pos3);
                plot(Z((years-1)*tf+1:years*tf),D.Shannon((years-1)*tf+1:years*tf),'.',Z((years-1)*tf+1:years*tf),G.Shannon((years-1)*tf+1:years*tf),'.',Z((years-1)*tf+1:years*tf),Totshannon((years-1)*tf+1:years*tf),'.');
                set(gca,'YTick', []);
                xlabel('Copepods abundance [µgC L^{-1}]');
                xlim([min(Z) max(Z)]);

            
        end

    end 
end  