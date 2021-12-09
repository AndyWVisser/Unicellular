
figure(21)
clf;
annualD = D.Btotot(end-365+1:end);
annualG = G.Btotot(end-365+1:end);
t = 1:365;
subplot(2,2,1)
annD = plot(t,annualD,'-r',1,annualD(1),'or');
subplot(2,2,2)
annG = plot(t,annualG,'-b',1,annualG(1),'ob');
subplot(2,2,3)
plotD = imagesc(log10(xo),y(:,1),squeeze(D.Bio(it,:,:)));
set(gca, 'YDir', 'normal','Clim',[0,10]);colormap(jet(32)); c=colorbar;

subplot(2,2,4)
plotG = imagesc(log10(xo),y(:,1),squeeze(G.Bio(it,:,:)));
set(gca, 'YDir', 'normal','Clim',[0,10]);colormap(jet(32)); c=colorbar;
for it = 1:365
    it
    set(plotD, 'CData', squeeze(D.Bio(it,:,:)));
    set(annD(2),'XData',it,'YData',annualD(it));
    set(plotG, 'CData', squeeze(G.Bio(it,:,:)));
    set(annG(2),'XData',it,'YData',annualG(it));
    drawnow
end