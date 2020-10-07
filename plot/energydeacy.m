clc;clear;

NX=703;
NZ=274;
L=20;
NT=8000;   

fp0=fopen('../cfspml2nd/output/energydecay.dat','rb');
energyviscocpml=fread(fp0,[NT],'float');
fclose(fp0);

fp1=fopen('../pml2nd/output/energydecay.dat','rb');
energyviscopml=fread(fp1,[NT],'float');
fclose(fp1);


figure(2);

dbviscocpml=todb(energyviscocpml,max(energyviscocpml));
dbviscopml=todb(energyviscopml,max(energyviscocpml));

it=1:NT;
plot(it,dbviscocpml,'b',it,dbviscopml,'r');
grid on;
legend('CPML','PML','Location','Best');

axis([0,NT,-160,20]);

set(gca,'xtick',[0:1000:8000],'ytick',[-160:20:20],'xticklabel',[0:1:8],'fontsize',12);

xlabel('Time (s)','fontsize',12);
ylabel('Energy Decay (dB)','fontsize',12);

set(gcf,'Position',[100 100 500 350]); 
set(gca,'Position',[.15 .15 .80 .80]); 

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.0 3.5]);
print ./Fig/energy.eps -depsc2 -r600;
