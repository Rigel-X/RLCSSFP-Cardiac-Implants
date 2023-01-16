%RECONSTRUCTION FOR 3D LC-SSFP

%set initial parameters
res=416;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\');
folder = 'C:\Users\Ricardo\Desktop\';
filepath = '';

%SPECIFY INPUT FILE 
% measfile = 'C:\Users\Ricardo\Desktop\meas_MID56_3D_radial_cine_LC_FID37504.dat'
% measfile = 'C:\Users\Ricardo\Desktop\3110_ES\meas_MID216_3D_LC_4_FID35793.dat'
% measfile = 'C:\Users\Ricardo\Desktop\meas_MID48_3DLC_acqs_FID37997.dat'
% measfile = 'C:\Users\Ricardo\Desktop\meas_MID21_3D_LC_4_FID37970.dat'
% measfile = 'C:\Users\Ricardo\Desktop\3Dphantoms\meas_MID192_3DTRY1_noint_FID38454.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID251_3D_R_192_15_16_nointnopc_nodyn_FID39973.dat'
measfile = 'C:\Users\Ricardo\Desktop\1205_YR\meas_MID85_3DTRY6_FID40560.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0812_bottlephantoms\meas_MID163_3DTRY_416_48x48_200bw_noPC_FID41271.dat'

%READ IN RAW DATA AND SORT, OUTPUTS FOUR MATRICES - ONE FOR EACH PASS OF
%SSFP
dispopt = 'on';
filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);
phases = sMDH.sLC.ushPhase;     %PHASES
coils= sMDH.sLC.ulChannelId;    %COILS
measure= sMDH.sLC.ushRepetition;%REPETITION
part = sMDH.sLC.ushPartition;   %SLICE
Np = sMDH.sLC.ushLine;          %PROJECTIONS
Nx = sMDH.ushSamplesInScan;     %SAMPLES
if mod(Np,2) == 0               %SPECIFY PROJECTION ANGLES GIVEN ODD OR EVEN Np
phival=[0:pi/Np:pi-pi/Np];
end
if mod(Np,2) ~=0
    phival=[0:2*pi/Np:2*pi-2*pi/Np];
end
nunder=1;% undersampling factor
phival= phival(1:nunder:end);
%SORT MATRICES INTO DESIRED FORMAT
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

%CLEAR VARIABLES TO SAVE MEMORY
clear rawdata loopcounters lc_names
clear first second third fourth
clear ksampsone ksampstwo ksampsthree ksampsfour

% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
k=[-.5:1/Nx:.5-1/Nx]*pi;
dcf1=abs(k).*cos(k);
dcf=repmat(dcf1,[Np 1]); 

%INITIALIZE MATRICES
ft_LC_0 = zeros(part, phases, coils, res, res); 
ft_LC_90 = zeros(part,phases, coils, res, res); 
ft_LC_180 = zeros(part,phases, coils, res, res); 
ft_LC_270 = zeros(part,phases, coils, res, res); 
% add = zeros(part,phases, coils, res, res); 
% sub = zeros(part,phases, coils, res, res); 
four = zeros(part,phases, coils, res, res); 

%create dcf-matrix (Ram-lak-filter) (same for all coildata)
% dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
% % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
% dcf2=dcf1/max(dcf1);
% dcf=repmat(dcf2,[Np 1]); %384*39 matrix
% figure;plot(dcf(1,:));
%end % if correct_all
%locate kspace samples in complex space (same for all coildata)

%CALCULATE RADIAL PROJECTIONS FOR GRIDDING
loc=zeros(measure,Nx,Np);
location =zeros(Nx,Np);
for m =1:measure
for j=[1:Nx] 
    for k=[1:Np] 
         if m == 1
             offset = 0;    %NO OFFSET IN FIRST PASS
         end
         if m > 1
             offset = (pi)/(measure*Np);    %CALCULATE INTERLEAVING OFFSET
         end
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
       
        pha=phival(k)+((m-1)*offset);
        loc(m,j,k)=mag*exp(1i*pha);
    end
end
end 

%PERFORM FFT IN Z DIRECTION (SLICE DIRECTION)
for ph = 1:phases, for cl = 1:coils, for i = 1:Nx, for j = 1:Np, 
        dataone(:,cl,ph,i,j) = ifftshift((fft(fftshift(ksampsnewone(:,cl,ph,i,j),1))),1); 
        datatwo(:,cl,ph,i,j) = ifftshift(fft(fftshift(ksampsnewtwo(:,cl,ph,i,j))));
        datathree(:,cl,ph,i,j) = ifftshift(fft(fftshift(ksampsnewthree(:,cl,ph,i,j))));
        datafour(:,cl,ph,i,j) = ifftshift(fft(fftshift(ksampsnewfour(:,cl,ph,i,j))));
end, end, end, end

%GRID EACH SLICE AND PERFORM 2D FFT WITHIN EACH SLICE, DEAPODIZATION
for pa=1:part,for cl = 1:coils,for ph = 1:phases,
            loctemp = squeeze(loc(1,:,:));
            datatemp = squeeze(dataone(pa,cl,ph,:,:));
            dcftemp = dcf.';
            
        datoneone(pa, cl,ph,:,:) = gridkb(loctemp(:),datatemp(:),dcftemp(:),res,kwidth,oversmpl);
        dattwoone(pa, cl,ph,:,:) = gridkb(loc(2,:,:),squeeze(squeeze(datatwo(pa,cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        datthreeone(pa, cl,ph,:,:) = gridkb(loc(3,:,:),squeeze(squeeze(datathree(pa,cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        datfourone(pa, cl,ph,:,:) = gridkb(loc(4,:,:),squeeze(squeeze(datafour(pa,cl,ph,:,:))),dcf',res,kwidth,oversmpl);
    
        datoneUP(pa, ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(datoneone(pa,cl,ph,:,:))))),oversmpl,res,kwidth);
        dattwoUP(pa, ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(dattwoone(pa, cl,ph,:,:))))),oversmpl,res,kwidth);
        datthreeUP(pa, ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(datthreeone(pa, cl,ph,:,:))))),oversmpl,res,kwidth);
        datfourUP(pa, ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(datfourone(pa, cl,ph,:,:))))),oversmpl,res,kwidth);

        end,end,end

%CLEAR NON IMPORTANT VARIABLES TO CONSERVE MEMORY
clear ksampsnewone ksampsnewtwo ksampsnewthree ksampsnewfour 
clear dataone datatwo datathree datafour datoneone dattwoone datthreeone datfourone
clear newraw

%PERFORM LINEAR COMBINATION (CAN BE DONE IN FOURIER SPACE)
%NOTE WHICH MATRIC IS SET TO FT_LC_9
four = zeros(part,phases, coils, res, res);
for pa=1:part
for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(pa,phase,coil,:,:) =(datoneUP(pa, phase, coil, :,:));
        ft_LC_90(pa,phase,coil,:,:) = (dattwoUP(pa, phase, coil, :,:));
        ft_LC_180(pa,phase,coil,:,:) = (datthreeUP(pa, phase, coil, :,:));
        ft_LC_270(pa,phase,coil,:,:) = (datfourUP(pa, phase, coil, :,:));
% 
% %         add(pa, phase, coil,:,:) = ft_LC_0(pa, phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(pa, phase, coil, :,:)));
% % % %         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
% %         sub(pa, phase, coil,:,:) = ft_LC_0(pa, phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(pa, phase, coil, :,:)));
% % %         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        four(pa, phase, coil,:,:) = (ft_LC_0(pa, phase, coil,:,:)) + (ft_LC_90(pa, phase, coil, :,:)) + (ft_LC_180(pa, phase, coil,:,:)) + (ft_LC_270(pa, phase, coil,:,:));
 end, end,end

% load ft_LC_0 
% load ft_LC_90 
% load ft_LC_180 
% load ft_LC_270

vidObj = VideoWriter('3D_four_LA2');
vidObj.FrameRate = 2;
open(vidObj);

% % phases = 1;
% phases = 1;
% 
% oot=1;
% phases = 1;

for pa = 1:part
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
% if oot == 1;
figure;finalIm(pa,phase,:,:)=mergeCoilImages(squeeze(squeeze((abs(four(pa,phase,:,:,:))))));
final_img(pa,phase,:,:) = (((abs(finalIm(pa,phase,:,:)))));
final(pa,phase,:,:) = fliplr(final_img(pa,phase,:,:));
    colormap(gray); imagesc(squeeze(final(pa,phase,:,:))); axis equal;  caxis([0,max(max(max(max(final))))]); title('FOUR STD 48p no PC');
set(gca,'Xdir','reverse')
% end
% % %  
% % % %     
% if oot==2;
% figure;finalIm(pa,phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(pa,phase,:,:,:))));
% final_img(pa,phase,:,:) = ((((finalIm(pa,phase,:,:)))));
% final(pa,phase,:,:) = fliplr(final_img(pa,phase,:,:));
%     colormap(gray); imagesc((squeeze((final(pa,phase,:,:))))); axis equal; caxis([0,max(max(max(max(final))))]); title('180 STD 64p');
%     set(gca,'Xdir','reverse')
% end
% if oot==4;
% figure;finalIm(pa,phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_0(pa,phase,:,:,:))));
% final_img(pa,phase,:,:) = ((((finalIm(pa,phase,:,:)))));
% final(pa,phase,:,:) = fliplr(final_img(pa,phase,:,:));
%     colormap(gray); imagesc((squeeze((final(pa,phase,:,:))))); axis equal; caxis([0,max(max(max(max((final)))))]); title('0 STD 64p');
%     set(gca,'Xdir','reverse')
% end
% if oot==4;
% figure;finalIm(pa,phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_90(pa,phase,:,:,:))));
% final_img(pa,phase,:,:) = ((((finalIm(pa,phase,:,:)))));
% final(pa,phase,:,:) = fliplr(final_img(pa,phase,:,:));
%     colormap(gray); imagesc((squeeze(flipud(final(pa,phase,:,:))))); axis equal; caxis([0,max(max(max(max(final))))]); title('90 STD 64p');
%     set(gca,'Xdir','reverse');
%     set(gca,'Ydir', 'reverse');
% end
% if oot==5;
% figure;finalIm(pa,phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_270(pa,phase,:,:,:))));
% final_img(pa,phase,:,:) = ((((finalIm(pa,phase,:,:)))));
% final(pa,phase,:,:) = fliplr(final_img(pa,phase,:,:));
%     colormap(gray); imagesc((squeeze((final(pa,phase,:,:))))); axis equal; caxis([0,max(max(max(max(final))))]); title('270 STD 64p');
%     set(gca,'Xdir','reverse')
% end

% 
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);

%% 
% % %% 
end
end
% voxel = [2 2 6];
% myimage=(permute(final,[3 4 1 2])); % create myimage\
% myimage = fliplr(myimage);
% myimage = flipud(myimage);
% 
% 
% if oot == 1;
%     imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms','3D_416_4_48_nopc'); %end %create a %%file name
% if oot == 2;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\3D_STD_INT_64\','3D_BOTTLES_STD_INT_64_180'); end
% if oot ==3;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\3D_STD_INT_64\','3D_BOTTLES_STD_INT_64_0'); end
% if oot == 4;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\3D_STD_INT_64\','3D_BOTTLES_STD_INT_64_90'); end
% if oot == 5;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\3D_STD_INT_64\','3D_BOTTLES_STD_INT_64_270'); end
% nii= make_nii(myimage, voxel);         
% save_nii(nii, imfilename);
% 
%   view_nii(nii);
%   clear nii
%   oot=oot+1;


close(vidObj);


vidObj = VideoWriter('3D_four_LV2_montage');
vidObj.FrameRate = 2;
open(vidObj);

montage = zeros(phases, 128,256);
for p = 1:phases, 
    montage(p, 1:64,1:64) = finalIm(3, p, 58:121, 38:101);
    montage(p, 1:64, 65:128) = finalIm(4, p, 58:121, 38:101);
    montage(p, 1:64,129:192) = finalIm(5, p, 58:121, 38:101);
    montage(p, 1:64,193:256) = finalIm(6, p, 58:121, 38:101);
    montage(p, 65:128,1:64) = finalIm(7, p, 58:121, 38:101);
	montage(p, 65:128,65:128) = finalIm(8, p, 58:121, 38:101);
    montage(p, 65:128,129:192) = finalIm(9, p, 58:121, 38:101);
    montage(p, 65:128,193:256) = finalIm(10, p, 58:121, 38:101);
    
    figure;colormap(gray); imagesc(fliplr(squeeze(abs(montage(p,:,:))))); axis equal;  caxis([0,max(max(max(max(montage))))]);
%     set(gca,'Xdir','reverse')
set(gca,'Visible','off')
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);