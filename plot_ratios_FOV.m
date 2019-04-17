clear all
close all

temp=load('15mindark_noiso_cit_FOV.mat');
dark_noiso=temp.FOV_ratios;
temp=load('15mindark_iso_cit_FOV.mat');
dark_iso=temp.FOV_ratios;
temp=load('15minfur_noiso_cit_FOV.mat');
fur_noiso=temp.FOV_ratios;
temp=load('15minfur_iso_cit_FOV.mat');
fur_iso=temp.FOV_ratios;
temp=load('15minlight_noiso_cit_FOV.mat');
light_noiso=temp.FOV_ratios;
temp=load('15minlight_iso_cit_FOV.mat');
light_iso=temp.FOV_ratios;

figure;
bar(1,mean(dark_noiso));
hold on
bar(2,mean(dark_iso));
bar(3,mean(fur_noiso));
bar(4,mean(fur_iso));
bar(5,mean(light_noiso));
bar(6,mean(light_iso));

errorbar(1,mean(dark_noiso),std(dark_noiso)./sqrt(length(dark_noiso)));
errorbar(2,mean(dark_iso),std(dark_iso)./sqrt(length(dark_iso)));
errorbar(3,mean(fur_noiso),std(fur_noiso)./sqrt(length(fur_noiso)));
errorbar(4,mean(fur_iso),std(fur_iso)./sqrt(length(fur_iso)));
errorbar(5,mean(light_noiso),std(light_noiso)./sqrt(length(light_noiso)));
errorbar(6,mean(light_iso),std(light_iso)./sqrt(length(light_iso)));

hold on
plot(1,dark_noiso,'o')
plot(2,dark_iso,'o')
plot(3,fur_noiso,'o')
plot(4,fur_iso,'o')
plot(5,light_noiso,'o')
plot(6,light_iso,'o')

fur_snr=mean(fur_iso)./mean(fur_noiso);
light_snr=mean(light_iso)./mean(light_noiso);

figure
cellmatrix=[dark_noiso dark_iso fur_noiso fur_iso light_noiso light_iso];
g1=[zeros(1,length(dark_noiso)) zeros(1,length(dark_iso)) ones(1,length(fur_noiso)) ones(1,length(fur_iso)) ones(1,length(light_noiso))*2 ones(1,length(light_iso))*2]; 
g2=[zeros(1,length(dark_noiso)) ones(1,length(dark_iso)) zeros(1,length(fur_noiso)) ones(1,length(fur_iso)) zeros(1,length(light_noiso)) ones(1,length(light_iso))]; 

[p,tbl,stats] = anovan(cellmatrix,{g1 g2},2);
multcompare(stats,'Dimension',[1 2])
