% plot response profile for bssfp with Partial dephasing 
% used for @ISMRM2022 poster presentation, jie xiang @yale mrrc

clear
close all
load('pulsatile_medium.mat')
d_phi = 2.5;
n = 4;
flow_max = 1.2; % cm/s

flow_max = max(find(flow_vect<flow_max));
color_param = linspace(0, 1, length(flow_vect(1:n:flow_max)))';
newDefaultColors = [color_param 1-color_param 1-color_param];

for i = 1:3
figure, set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
plot(d_phi:d_phi:360, (spectral_profs(1:n:flow_max, :, i))','LineWidth',2)
xlabel('Precession per TR (^o)'), ylabel('|M_{xy}| (a.u.)')
axis([0 360 0 18])
end

spin_replacement_rate = 100 * flow_vect * 3.5/6; % 3.5 ms TR, 6 mm FWHM slice
spin_replacement_rate = 100 * flow_vect(1:n:flow_max); % cm/s
figure, imagesc(round(spin_replacement_rate)), colormap(newDefaultColors), cbr_h = colorbar;
set(cbr_h, 'YTick', round(spin_replacement_rate(1:2:end)))


%%
newDefaultColors =[
1 0.85 0.1;
0.9290 0.6940 0.1250;
0.8990 0.5040 0.1120;
0.8500 0.3250 0.0980;
0.7400 0.2050 0.1480;
0.6350 0.0780 0.2040
0.5200 0      0.2500]

%%
color_param = linspace(0.3, 1, 7)';
% newDefaultColors = [color_param 1-color_param 1-color_param];
figure, 
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
plot(d_phi:d_phi:360, (spectral_profs(1, :, 1))','LineWidth',2),hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 2))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 3))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 4))'./(max(spectral_profs(1, :, 4))),'LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 5))'./(max(spectral_profs(1, :, 5))),'LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 6))'./(max(spectral_profs(1, :, 6))),'LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(1, :, 7))'./(max(spectral_profs(1, :, 7))),'LineWidth',2)
xlabel('Precession per TR (^o)'), ylabel('|M_{xy}| (a.u.)')
axis([0 360 0 2])
ylim([0.1,1.1])

%%
color_param = linspace(0, 1, 7)';
% newDefaultColors = [color_param 1-color_param 1-color_param];
figure, 
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
plot(d_phi:d_phi:360, (spectral_profs(end, :, 1))','LineWidth',2),hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 2))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 3))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 4))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 5))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 6))','LineWidth',2), hold on
plot(d_phi:d_phi:360, (spectral_profs(end, :, 7))'./(max(spectral_profs(end, :, 7))),'LineWidth',2)
xlabel('Precession per TR (^o)'), ylabel('|M_{xy}| (a.u.)')
axis([0 360 0 18])
ylim([0,12])
legend('0','1/12','2/12','3/12','4/12','5/12','6/12')

%%
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
% ylim([0,10])

lgd = legend('0','1/12','2/12','3/12','4/12','5/12');
lgd.Color = 'black';
lgd.TextColor = 'white';

%%
grid on
lgd = legend('0','1/12','2/12','3/12','4/12','5/12','6/12');
