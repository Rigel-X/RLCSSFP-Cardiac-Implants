% gridding based on hargreaves spiral gridding code.


folder = 'C:\Users\Dana Peters\Documents\';
filepath = '';

%PHANTOMFORCOMPARISON
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID47_LC_C12_0_125_FID30409.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID48_LC_C12_90_125_FID30410.dat'
measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID49_LC_C12_180_125_FID30411.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID50_LC_C12_270_125_FID30412.dat'

%HUMAN SCAN FILES
% measfile = 'C:\Users\Ricardo\Desktop\0906_YR_twix\meas_MID102_0830_bSSFP_C_180_FID26224.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0906_YR_twix\meas_MID104_0830_bSSFP_C_90_FID26226.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0906_YR_twix\meas_MID103_0830_bSSFP_C_0_FID26225.dat'
%  measfile = 'C:\Users\Ricardo\Desktop\0906_YR_twix\meas_MID105_0830_bSSFP_C_270_FID26227.dat'
%HUMAN SCAN FILES 2.0
% measfile = 'C:\Users\Ricardo\Desktop\0914_YR2twix\meas_MID113_LC_Cart_0_FID28357.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0914_YR2twix\meas_MID114_LC_Cart_90_FID28358.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0914_YR2twix\meas_MID112_LC_Cart_180_FID28356.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0914_YR2twix\meas_MID115_LC_Cart_270_FID28359.dat'
%PHANTOM FILES
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID23_0830_bSSFP_C_0_2_FID27775.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID22_0830_bSSFP_C_90_FID27774.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID19_0830_bSSFP_C_180_FID27771.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID24_0830_bSSFP_C_270_FID27776.dat'
%PHANTOM FILES - 50 Offset
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID34_0830_bSSFP_C_0_plus50_FID27786.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID35_0830_bSSFP_C_90_plus50_FID27787.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID33_0830_bSSFP_C_180_plus50_FID27785.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID36_0830_bSSFP_C_270_plus50_FID27788.dat'
%PHANTOM FILES - 100 Offset
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID37_0830_bSSFP_C_0_plus100_FID27789.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID39_0830_bSSFP_C_90_plus100_FID27791.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID38_0830_bSSFP_C_180_plus100_FID27790.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0911_PHAN\meas_MID40_0830_bSSFP_C_270_plus100_FID27792.dat'
%RG
% measfile = 'C:\Users\Ricardo\Desktop\0921_RG\meas_MID145_LC_C_0_FID29438.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0921_RG\meas_MID146_LC_C_90_FID29439.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0921_RG\meas_MID144_LC_C_180_FID29437.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0921_RG\meas_MID147_LC_C_270_FID29440.dat'

% measfile = 'C:\Users\Ricardo\Desktop\0926_phan\meas_MID90_LC_C12_180_75_FID29964.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0926_phan\meas_MID89_LC_C12_0_75_FID29963.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0926_phan\meas_MID88_LC_C12_90_FID29962.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0926_phan\meas_MID91_LC_C12_270_75_FID29965.dat'

% measfile = 'C:\Users\Ricardo\Desktop\0926_phan\meas_MID78_LC_C12_180_75FID29964.dat'
dispopt = 'on';
[rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17phasenew(measfile, dispopt);
phases = sMDH.sLC.ushPhase;         %phases
nc=sMDH.sLC.ulChannelId;            %coils
npe = sMDH.sLC.ushLine;             %lines
ns = sMDH.ushSamplesInScan;         %samples
centre = ns/2;
%newraw (coil, phase, line, sample), (32, 6, 126, 256)
ft_C_oe = zeros(nc, phases, npe, ns);
% ft_LC = zeros(nc,phases,npe,ns);
% add = zeros(nc, phases, npe, ns);
% sub = zeros(nc, phases, npe, ns);
% four = zeros(nc, phases, npe, ns);
% sum_add = zeros(nc, phases, npe, ns);
% sum_sub = zeros(nc, phases, npe, ns);
% sum_four = zeros(nc, phases, npe, ns);


% vidObj = VideoWriter('ADD');
% open(vidObj);

for coil = 1:nc
    for phase = 1:phases
        ft_C_oe(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(newraw(coil, phase, :, :)))));
%         add(coil,phase,:,:) = ft_C_ze(coil, phase,:,:)+((sqrt(-1))*(ft_C_oe(coil, phase, :,:)));
%         sub(coil,phase,:,:) = ft_C_ze(coil, phase,:,:)-((sqrt(-1))*(ft_C_oe(coil, phase, :,:)));
%         four(coil,phase,:,:) = (ft_C_ze(coil, phase,:,:)) + (ft_C_oe(coil, phase, :,:)) + (ft_C_ts(coil, phase,:,:)) + (ft_C_nn(coil, phase, :,:));
    end
end
% 
% phase =1;
% for coil = 1:nc

% %     figure;imagesc(squeeze(angle(add(coil,phase,:,:)))); axis equal; title(sprintf( "%s%s%d", 'ADD Magnitude  ', 'Coil = ', nc));
% %     figure;imagesc(squeeze(abs(add(coil,phase,:,:)))); axis equal; title(sprintf( "%s%s%d", 'ADD Phase  ', 'Coil = ', nc));
% % %         
% end
phases = 1;
for phase = 1:phases
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(add(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD Phase   ');
% % % % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(add(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD   ');
% % % %     
% % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(sub(:,phase,:,:))));
% % % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('SUBTRACT Phase  ');
% % % % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(sub(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('SUBTRACT   ');
% % % %     
% % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(four(:,phase,:,:))));
% % % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('FOUR Phase ) ');
% % % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('FOUR   ');
% %     
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_C_oe(:,phase,:,:))));
    colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('OE  ');     

% % % 

%         currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
end

% close(vidObj);
