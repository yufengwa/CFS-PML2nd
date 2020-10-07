clc;clear;

NX=703;
NZ=274;
L=20;
NT=8000;


Filename = 'snapshots.gif';
set(0,'defaultfigurecolor','w');


fp1=fopen('../pml2nd/output/visco_1_1_snapshots1_800.dat','rb');
snap1=fread(fp1,[NZ NX],'float');
fclose(fp1);

fp2=fopen('../pml2nd/output/visco_1_1_snapshots1_1600.dat','rb');
snap11=fread(fp2,[NZ NX],'float');
fclose(fp2);


fp3=fopen('../pml2nd/output/visco_1_1_snapshots1_2800.dat','rb');
snap111=fread(fp3,[NZ NX],'float');
fclose(fp3);


fp4=fopen('../pml2nd/output/visco_1_1_snapshots1_4800.dat','rb');
snap1111=fread(fp4,[NZ NX],'float');
fclose(fp4);

fp5=fopen('../pml2nd/output/visco_1_1_snapshots1_6800.dat','rb');
snap11111=fread(fp5,[NZ NX],'float');
fclose(fp5);
figure;

    imagesc(1e8*snap1);
    hold on;
    
    %plot(270,28,'bp',62,92,'s',212,172,'s',462,212,'s',652,42,'s','MarkerSize',8,'MarkerFaceColor',[1,0.8,0]);
    plot(270,28,'bp','MarkerSize',12,'MarkerFaceColor',[1,0.8,0]);
    hold on;
    
    line([L,L],[0,NZ],'Color',[1,0,0]);
    line([NX-L,NX-L],[0,NZ],'Color',[1,0,0]);
    line([0,NX],[L,L],'Color',[1,0,0]);
    line([0,NX],[NZ-L,NZ-L],'Color',[1,0,0]);
    hold on;
    
    axis ([0,NX,0,NZ]);
    set(gca,'xtick',[0:100:700],'ytick',[0:100:300],'xticklabel',[0:1:7],'yticklabel',[0:1:3],'fontsize',10);
    colormap(gray);
    mina=min(caxis);
    maxa=max(caxis);
    caxis([0.1*mina,0.1*maxa]);
    colorbar('EastOutside');
    axis image;
    xlabel('Distance (km)');
    ylabel('Depth (km)');

    set(gcf,'Position',[100 100 640 280]); 
    set(gca,'Position',[.07 .08 .78 .96]); 

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.4 2.8]);
    print ./Fig/pmlsnapmar1.eps -depsc2 -r600;

    
    
figure;

    imagesc(1e8*snap11); 
    hold on;
    
    
    %plot(270,28,'bp',62,92,'s',212,172,'s',462,212,'s',652,42,'s','MarkerSize',8,'MarkerFaceColor',[1,0.8,0]);
    plot(270,28,'bp','MarkerSize',12,'MarkerFaceColor',[1,0.8,0]);
    hold on;
    
    line([L,L],[0,NZ],'Color',[1,0,0]);
    line([NX-L,NX-L],[0,NZ],'Color',[1,0,0]);
    line([0,NX],[L,L],'Color',[1,0,0]);
    line([0,NX],[NZ-L,NZ-L],'Color',[1,0,0]);
    hold on;
    
    axis ([0,NX,0,NZ]);
    set(gca,'xtick',[0:100:700],'ytick',[0:100:300],'xticklabel',[0:1:7],'yticklabel',[0:1:3],'fontsize',10);
    colormap(gray);
    mina=min(caxis);
    maxa=max(caxis);
    caxis([0.1*mina,0.1*maxa]);
    colorbar('EastOutside');
    axis image;
    xlabel('Distance (km)');
    ylabel('Depth (km)');

    set(gcf,'Position',[100 100 640 280]); 
    set(gca,'Position',[.07 .08 .78 .96]); 

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.4 2.8]);
    print ./Fig/pmlsnapmar2.eps -depsc2 -r600;
    
    
figure;

    imagesc(1e8*snap111); 
    hold on;
    
    %plot(270,28,'bp',62,92,'s',212,172,'s',462,212,'s',652,42,'s','MarkerSize',8,'MarkerFaceColor',[1,0.8,0]);
    plot(270,28,'bp','MarkerSize',12,'MarkerFaceColor',[1,0.8,0]);
    hold on;
    
    line([L,L],[0,NZ],'Color',[1,0,0]);
    line([NX-L,NX-L],[0,NZ],'Color',[1,0,0]);
    line([0,NX],[L,L],'Color',[1,0,0]);
    line([0,NX],[NZ-L,NZ-L],'Color',[1,0,0]);
    hold on;
    
    axis ([0,NX,0,NZ]);
    set(gca,'xtick',[0:100:700],'ytick',[0:100:300],'xticklabel',[0:1:7],'yticklabel',[0:1:3],'fontsize',10);
    colormap(gray);
    mina=min(caxis);
    maxa=max(caxis);
    caxis([mina,0.6*maxa]);
    colorbar('EastOutside');
    axis image;
    xlabel('Distance (km)');
    ylabel('Depth (km)');

    set(gcf,'Position',[100 100 640 280]); 
    set(gca,'Position',[.07 .08 .78 .96]); 

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.4 2.8]);
    print ./Fig/pmlsnapmar3.eps -depsc2 -r600;

    
    
figure;

    imagesc(1e8*snap1111); 
    hold on;
    
    %plot(270,28,'bp',62,92,'s',212,172,'s',462,212,'s',652,42,'s','MarkerSize',8,'MarkerFaceColor',[1,0.8,0]);
    plot(270,28,'bp','MarkerSize',12,'MarkerFaceColor',[1,0.8,0]);
    hold on;
    
    line([L,L],[0,NZ],'Color',[1,0,0]);
    line([NX-L,NX-L],[0,NZ],'Color',[1,0,0]);
    line([0,NX],[L,L],'Color',[1,0,0]);
    line([0,NX],[NZ-L,NZ-L],'Color',[1,0,0]);
    hold on;
    
    axis ([0,NX,0,NZ]);
    set(gca,'xtick',[0:100:700],'ytick',[0:100:300],'xticklabel',[0:1:7],'yticklabel',[0:1:3],'fontsize',10);
    colormap(gray);
    mina=min(caxis);
    maxa=max(caxis);
    caxis([1.2*mina,1.2*maxa]);
    colorbar('EastOutside');
    axis image;
    xlabel('Distance (km)');
    ylabel('Depth (km)');

    set(gcf,'Position',[100 100 640 280]); 
    set(gca,'Position',[.07 .08 .78 .96]); 

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.4 2.8]);
    print ./Fig/pmlsnapmar4.eps -depsc2 -r600;


figure;

    imagesc(1e8*snap11111); 
    hold on;
    
    %plot(270,28,'bp',62,92,'s',212,172,'s',462,212,'s',652,42,'s','MarkerSize',8,'MarkerFaceColor',[1,0.8,0]);
    plot(270,28,'bp','MarkerSize',12,'MarkerFaceColor',[1,0.8,0]);
    hold on;
    
    line([L,L],[0,NZ],'Color',[1,0,0]);
    line([NX-L,NX-L],[0,NZ],'Color',[1,0,0]);
    line([0,NX],[L,L],'Color',[1,0,0]);
    line([0,NX],[NZ-L,NZ-L],'Color',[1,0,0]);
    hold on;
    
    axis ([0,NX,0,NZ]);
    set(gca,'xtick',[0:100:700],'ytick',[0:100:300],'xticklabel',[0:1:7],'yticklabel',[0:1:3],'fontsize',10);
    colormap(gray);
    mina=min(caxis);
    maxa=max(caxis);
    caxis([1.2*mina,1.2*maxa]);
    colorbar('EastOutside');
    axis image;
    xlabel('Distance (km)');
    ylabel('Depth (km)');

    set(gcf,'Position',[100 100 640 280]); 
    set(gca,'Position',[.07 .08 .78 .96]); 

    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.4 2.8]);
    print ./Fig/pmlsnapmar5.eps -depsc2 -r600;    