% gridding based on hargreaves spiral gridding code.


folder = 'C:\Users\Dana Peters\Documents\';
filepath = '';
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID68_LC_offset_R30_without_FID30891.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID70_LC_offset_R60_without_FID30893.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID84_LC_offset_R30_without_125_FID30907.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID88_LC_offset_R60_without_125_FID30911.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID62_LC_C_75_60fp_FID30423.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID66_LC_C_75_60fp_30seg_FID30427.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID56_LC_C_75_FID30417.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID45_LC_C_125_FID30407.dat'
measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID59_LC_C_FID31667.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID77_LC_C_SHAX_FID31684.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID78_LC_C_SHAX_TR3_8_FID31685.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID79_LC_C_SHAX_TR5_7_FID31686.dat'
dispopt = 'on';
[first, second, third, fourth, rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);

phases = sMDH.sLC.ushPhase;         %phases
nc=sMDH.sLC.ulChannelId;            %coils
npe = sMDH.sLC.ushLine;             %lines
ns = sMDH.ushSamplesInScan;         %samples
centre = ns/2;
%newraw (coil, phase, line, sample), (32, 6, 126, 256)
ft_C_ts = zeros(nc, phases, npe, ns);
ft_LC = zeros(nc,phases,npe,ns);
add = zeros(nc, phases, npe, ns);
sub = zeros(nc, phases, npe, ns);
four = zeros(nc, phases, npe, ns);
sum_add = zeros(nc, phases, npe, ns);
sum_sub = zeros(nc, phases, npe, ns);
sum_four = zeros(nc, phases, npe, ns);


vidObj = VideoWriter('LC_LA_180');
open(vidObj);

for coil = 1:nc
    for phase = 1:phases
        ft_LC_0(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(first(coil, phase, :, :)))));
        ft_LC_90(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(second(coil, phase, :, :)))));
        ft_LC_180(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(third(coil, phase, :, :)))));
        ft_LC_270(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(fourth(coil, phase, :, :)))));
        add(coil,phase,:,:) = ft_LC_0(coil, phase,:,:)+((sqrt(-1))*(ft_LC_180(coil, phase, :,:)));
        sub(coil,phase,:,:) = ft_LC_0(coil, phase,:,:)-((sqrt(-1))*(ft_LC_180(coil, phase, :,:)));
        four(coil,phase,:,:) = (ft_LC_0(coil, phase,:,:)) + (ft_LC_180(coil, phase, :,:)) + (ft_LC_270(coil, phase,:,:)) + (ft_LC_90(coil, phase, :,:));
    end
end
% 
% phase =1;
% for coil = 1:nc

% %     figure;imagesc(squeeze(angle(add(coil,phase,:,:)))); axis equal; title(sprintf( "%s%s%d", 'ADD Magnitude  ', 'Coil = ', nc));
% %     figure;imagesc(squeeze(abs(add(coil,phase,:,:)))); axis equal; title(sprintf( "%s%s%d", 'ADD Phase  ', 'Coil = ', nc));
% % %         
% end
% phases = 1;
for phase = 1:phases
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(add(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD Phase   ');
% % % % % % 
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
% % % % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(:,phase,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('FOUR   ');
% % %     
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(:,phase,:,:))));
    colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('180  ');     

% % % 

        currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end

close(vidObj);
