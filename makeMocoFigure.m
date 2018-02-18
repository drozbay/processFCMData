clear all; close all;

load dataForMOCOFig.mat;

mocoData_11 = importMOCO(path11,infoFCM_11);
mocoData_13 = importMOCO(path13,infoFCM_13);
mocoData_15 = importMOCO(path15,infoFCM_15);
mocoData_17 = importMOCO(path17,infoFCM_17);

figure(1); clf;
subplot(2,1,1);hold on;
plot(mocoData_11.t_s,mocoData_11.mag_um);
plot(mocoData_13.t_s,mocoData_13.mag_um);
plot(mocoData_15.t_s,mocoData_15.mag_um);
plot(mocoData_17.t_s,mocoData_17.mag_um);
subplot(2,1,2); hold on;
gradN = 1;
plot(mocoData_11.t_s,gradient(mocoData_11.mag_um,gradN));
plot(mocoData_13.t_s,gradient(mocoData_13.mag_um,gradN));
plot(mocoData_15.t_s,gradient(mocoData_15.mag_um,gradN));
plot(mocoData_17.t_s,gradient(mocoData_17.mag_um,gradN));

mocoMag_all = [mocoData_11.mag_um; mocoData_13.mag_um; mocoData_15.mag_um; mocoData_17.mag_um;];
mocoMag_allDiff = [gradient(mocoData_11.mag_um,gradN); gradient(mocoData_13.mag_um,gradN); gradient(mocoData_15.mag_um,gradN); gradient(mocoData_17.mag_um,gradN);];

totalDuration = mocoData_11.t_s(end)+mocoData_13.t_s(end)+mocoData_15.t_s(end)+mocoData_17.t_s(end);

figure(2);
histEdges = (0:1:10);
subplot(1,2,1);
h1 = histogram(mocoMag_all,histEdges);
subplot(1,2,2);
h2 = histogram(mocoMag_allDiff,histEdges/4);

figure(3);
histEdges = (0:1:10);
h3 = histogram(mocoMag_all,histEdges);
h3Counts = histcounts(mocoMag_all,histEdges);
h3.FaceColor = [0.2 0.2 0.2];
title({'Cumulative motion artifacts',...
    sprintf('Total imaging duration: %d:%d minutes',...
    floor(totalDuration/60),ceil(rem(totalDuration,60)))});
xlabel('Lateral motion artifact (\mum)');
ylabel('Counts');
set(gca,'Fontsize',8);
set(gcf,'units','inches');
curPos = get(gcf,'Position');
set(gcf,'Position',[curPos(1),curPos(2),3,3]);
t1 = text(histEdges(end)-1,max(h3Counts(:))*0.9,...
    sprintf('Mean = %0.1f \\mum\nStdDev = %0.1f \\mum',...
    mean(mocoMag_all),std(mocoMag_all,1)),...
    'HorizontalAlignment','right','VerticalAlignment','Top');
