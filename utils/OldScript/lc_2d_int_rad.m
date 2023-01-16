% gridding based on hargreaves spiral gridding code.

%set initial parameters
res=416;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\');
folder = 'C:\Users\Ricardo\Desktop\';
filepath = '';

measfile='C:\Users\Ricardo\Desktop\0812_bottlephantoms\meas_MID136_LC_R_GX_COMP_INT_48_FID41244.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID52_LC_R_018090270_SL1_FID36499.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID53_LC_R_027090180_SL1_FID36500.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID62_LC_R_0090180270_SL2_FID36509.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID63_LC_R_018090270_SL2_FID36510.dat'


dispopt = 'on';
filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);
phases = sMDH.sLC.ushPhase;
coils= sMDH.sLC.ulChannelId; %sMDH.sLC.ushPartition; sMDH.ushUsedChannels;
measure= sMDH.sLC.ushRepetition;
part = sMDH.sLC.ushPartition;
Np = sMDH.sLC.ushLine;
Nx = sMDH.ushSamplesInScan;
phival=[0:pi/Np:pi-pi/Np];
nunder=1;% undersampling factor
phival= phival(1:nunder:end);
ksampsone = first;
ksampsnewone=permute(ksampsone,[1 2 3 5 4]);
ksampsnewone= ksampsnewone(:,:,:,:,1:nunder:end);
ksampstwo = second;
ksampsnewtwo=permute(ksampstwo,[1 2 3 5 4]);
ksampsnewtwo= ksampsnewtwo(:,:,:,:,1:nunder:end);
ksampsthree = third;
ksampsnewthree=permute(ksampsthree,[1 2 3 5 4]);
ksampsnewthree= ksampsnewthree(:,:,:,:,1:nunder:end);
ksampsfour = fourth;
ksampsnewfour=permute(ksampsfour,[1 2 3 5 4]);
ksampsnewfour= ksampsnewfour(:,:,:,:,1:nunder:end);

% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
k=[-.5:1/Nx:.5-1/Nx];
dcf1=abs(k).*sinc(k);
dcf=repmat(dcf1,[Np 1]); 

ft_LC_0 = zeros(phases, coils, res, res); 
ft_LC_90 = zeros(phases, coils, res, res); 
ft_LC_180 = zeros(phases, coils, res, res); 
ft_LC_270 = zeros(phases, coils, res, res); 
add = zeros(phases, coils, res, res); 
sub = zeros(phases, coils, res, res); 
four = zeros(phases, coils, res, res); 

%create dcf-matrix (Ram-lak-filter) (same for all coildata)
% dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
% % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
% dcf2=dcf1/max(dcf1);
% dcf=repmat(dcf2,[Np 1]); %384*39 matrix
% figure;plot(dcf(1,:));

%end % if correct_all
%locate kspace samples in complex space (same for all coildata)
loc=zeros(measure,Nx,Np);
location =zeros(Nx,Np);
% C = lines(measure);

for m =1:measure
for j=[1:Nx] 
    for k=[1:Np] 
         if m == 1
             offset = 0;
         end
         if m > 1
             offset = (pi)/(measure*Np);
         end
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
       
        pha=phival(k)+((m-1)*offset);
        loc(m,j,k)=mag*exp(1i*pha);

    end
end
% C = {'c' 'b' 'g' 'm'};
% location = squeeze(loc(m,:,:));
% plot(real(location),imag(location), 'color',C(m));
% hold one
end
% h1 = plot(real(squeeze(loc(1,:,1))),imag(squeeze(loc(1,:,1))), 'color','c');hold on 
% h2 = plot(real(squeeze(loc(2,:,1))),imag(squeeze(loc(2,:,1))), 'color','b');
% h3 = plot(real(squeeze(loc(3,:,1))),imag(squeeze(loc(3,:,1))), 'color','g');
% h4 = plot(real(squeeze(loc(4,:,1))),imag(squeeze(loc(4,:,1))), 'color','m');
% 
% legend
% legend([h1 h2 h3 h4],{'First Pass','Second Pass','Third Pass','Fourth Pass'}, 'Box', 'off')
% h5= plot(real(squeeze(loc(1,:,2:Np))),imag(squeeze(loc(1,:,2:Np))), 'color','c', 'HandleVisibility','off'); 
% h6= plot(real(squeeze(loc(2,:,2:Np))),imag(squeeze(loc(2,:,2:Np))), 'color','b', 'HandleVisibility','off');
% h7= plot(real(squeeze(loc(3,:,2:Np))),imag(squeeze(loc(3,:,2:Np))), 'color','g', 'HandleVisibility','off');
% h8= plot(real(squeeze(loc(4,:,2:Np))),imag(squeeze(loc(4,:,2:Np))), 'color','m', 'HandleVisibility','off');
% set(gca,'Visible','off')


for pa=1:part
for cl = 1:coils
    for ph = 1:phases
        size(squeeze(ksampsnewone(pa,cl,ph,:,:)))
        datoneone(ph,cl,:,:) = gridkb(loc(1,:,:),squeeze(ksampsnewone(pa,cl,ph,:,:)),dcf',res,kwidth,oversmpl);
        dattwoone(ph,cl,:,:) = gridkb(loc(2,:,:),squeeze(ksampsnewtwo(pa,cl,ph,:,:)),dcf',res,kwidth,oversmpl);
        datthreeone(ph,cl,:,:) = gridkb(loc(3,:,:),squeeze(ksampsnewthree(pa,cl,ph,:,:)),dcf',res,kwidth,oversmpl);
        datfourone(ph,cl,:,:) = gridkb(loc(4,:,:),squeeze(ksampsnewfour(pa,cl,ph,:,:)),dcf',res,kwidth,oversmpl);
        imone(ph,cl,:,:) = datoneone(ph,cl,:,:) + dattwoone(ph,cl,:,:) + datthreeone(ph,cl,:,:)+ datfourone(ph,cl,:,:);
%         imone(ph,cl,:,:) = (fft2(squeeze(datone(ph,cl,:,:))));
%         imtwo(ph,cl,:,:) = (fft2(squeeze(dattwo(ph,cl,:,:))));
%         imthree(ph,cl,:,:) = (fft2(squeeze(datthree(ph,cl,:,:))));
%         imfour(ph,cl,:,:) = (fft2(squeeze(datfour(ph,cl,:,:))));
    end
%     finalIm(ph,:,:)=mergeCoilImages(squeeze(abs(im(ph,:,:,:))));
%     colormap(gray); imagesc(fftshift(squeeze(finalIm(ph,:,:))));
end 
end


for coil = 1:coils
    for phase = 1:phases
        four(phase,coil,:,:) = deapod(fftshift(fft2(fftshift(squeeze(imone(ph,cl,:,:))))),1,res,kwidth);
%         ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
%         ft_LC_90(phase, coil,:,:) = (imtwo(phase, coil, :, :));
%         ft_LC_180(phase, coil,:,:) = (imthree(phase, coil, :, :));
%         ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
%         add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
%         four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
    end
end

% 
% vidObj = VideoWriter('TRIAL');
% vidObj.FrameRate = 2;
% open(vidObj);
% 
% % % phases = 1;
% % phases = 1;
% 
for phase = 1:phases
% % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(add(phase,:,:,:))));
% % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD Phase   ');
% % % % % % % % % 
% %  figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(add(phase,:,:,:))));
% % final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
% % final(phase,:,:) = flip(final_img(phase,:,:),3);
% %     colormap(gray); imagesc(squeeze(final(phase,:,:))))); axis equal; title('ADD180270 - 090180270 - with offset  ');
% % % % % % % 
% % % % % % % % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(sub(phase,:,:,:))));
% % % % % % % % % % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('SUBTRACT Phase  ');
% % % % % % % 
% % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(sub(phase,:,:,:))));
% % final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
% % final(phase,:,:) = flip(final_img(phase,:,:),3);
% %     colormap(gray); imagesc(squeeze(final(phase,:,:))))); axis equal; title('SUBTRACT - 018090270 - with offset  ');
% % % % % % % % % % %     
% % % % % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(four(phase,:,:,:))));
% % % % % % %     colormap(gray); imagesc(fliplr(((squeeze(finalIm(phase,:,:))))); axis equal; title('FOUR Phase ) ');
% 
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = flip(final_img(phase,:,:),3);
    colormap(gray); imagesc(squeeze(final(phase,:,:))); axis equal; caxis([0,max(max(max(final)))]); title('FOUR COMP INT 48p');
% %  
% % %     
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(phase,:,:,:))));
% final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
% final(phase,:,:) = flip(final_img(phase,:,:),3);
%     colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; title('180');  
% %     set(gca,'Xdir','reverse')
% %     
% %     figure; colormap(gray);imagesc(fftshift(squeeze(ft_LC_180(phase,coil,:,:)))); axis equal;title('180');
%     
%         currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
% 

end
%% 
% % %% 

% voxel = [2 2 8];
% myimage=(permute(final,[3 2 1])); % create myimage\
% % myimage = flip(myimage, 2);
% 
% imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\1411_MM\','MM_LC_STD_four'); %create a %%file name
% 
% nii= make_nii(myimage, voxel);         
% save_nii(nii, imfilename);

%   view_nii(nii);
%   clear nii
% end
% close(vidObj);
% 
% 
% vidObj = VideoWriter('FOUR');
% open(vidObj);
% phases = 19;
% for phase=1:phases
%     figure
%     subplot(2,2,1)
%     four = zeros(phases,6, 160,160);
%     load four_IO
%     finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
%     final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
%     final(phase,:,:) = flip(final_img(phase,:,:),3);
%     colormap(gray); imagesc((squeeze((final(phase,:,:))))); caxis([0,max(max(max(final)))]);axis equal; title('IO'); 
%     
%     subplot(2,2,2)
%     four = zeros(phases,6, 160,160);
%     load four_DP
%     four(16,:,:,:) = four(15,:,:,:);four(17,:,:,:) = four(15,:,:,:);four(18,:,:,:) = four(15,:,:,:);four(19,:,:,:) = four(15,:,:,:);
%     finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
%     final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
%     final(phase,:,:) = flip(final_img(phase,:,:),3);
%     colormap(gray); imagesc((squeeze((final(phase,:,:))))); caxis([0,max(max(max(final)))]);axis equal; title('DP');
%     
%     subplot(2,2,3)
%     four = zeros(phases,6, 160,160);
%     load four_JA
%     four(17,:,:,:) = four(16,:,:,:);four(18,:,:,:) = four(17,:,:,:); four(19,:,:,:) = four(17,:,:,:);
%     finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
%     final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
%     final(phase,:,:) = flip(final_img(phase,:,:),3);
%     colormap(gray); imagesc((squeeze((final(phase,:,:))))); caxis([0,max(max(max(final)))]);axis equal; title('JA');
%     
%     subplot(2,2,4)
%     four = zeros(phases,6, 160,160);
%     load four_KR
%     four(19,:,:,:) = four(18,:,:,:);
%     finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
%     final_img(phase,:,:) = ((fftshift((finalIm(phase,:,:)))));
%     final(phase,:,:) = flip(final_img(phase,:,:),3);
%     colormap(gray); imagesc((squeeze((final(phase,:,:))))); caxis([0,max(max(max(final)))]);axis equal; title('KR');
%     
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
% 
% % voxel = [2 2 8];
% % myimage=(permute(final,[3 2 1])); % create myimage\
% % % myimage = flip(myimage, 2);
% % 
% % imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\gridding_orange_3T'); %create a %%file name
% % 
% % nii= make_nii(h, voxel);         
% % save_nii(nii, imfilename);
% % 
% %   view_nii(nii);
% %   clear nii
% end
% % 
% close(vidObj);