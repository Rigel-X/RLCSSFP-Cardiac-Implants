

%SET INITIAL PARAMTERS 
res=416;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\');
folder = 'C:\Users\Ricardo\Desktop\';
filepath = '';

% measfile = 'C:\Users\Ricardo\Desktop\1022_phans\timtriob_20181025_064338988\timtriob_20181022_164901069\meas_MID91_LC_R_018090270_80p_FID34336.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID256_LC_R_0090180270_FID34977.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID257_LC_R_018090270_FID34978.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID258_LC_R_027090180_FID34979.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_045134283\meas_MID210_3D_LC_FID34933.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_045134283\meas_MID212_3D_LC_4_FID34935.dat'

% %perps
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID278_LC_R_0090180270_FID34999.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID279_LC_R_018090270_P_FID35000.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID281_LC_R_027090180_Pcor_FID35002.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1011_phans\timtrio_20181011_164901069\meas_MID52_LC_R_0090180270_GX75_FID34297.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1011_phans\timtrio_20181011_164901069\meas_MID53_LC_R_018090270_GX75_FID34298.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1011_phans\timtrio_20181011_164901069\meas_MID54_LC_R_027090180_GX75_FID34299.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\timtriob_20181018_161230220\meas_MID239_LC_fresh_R_32_sl_FID33651.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\timtriob_20181018_161230220\meas_MID241_LC_fresh_R_48_sl_FID33653.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\timtriob_20181018_161230220\meas_MID243_LC_fresh_R_64_sl_FID33655.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID316_LC_fresh_R_64_sl_FID33728.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID318_LC_fresh_R_48_sl_FID33730.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID320_LC_fresh_R_32_sl_FID33732.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID325_LC_fresh_R_64_sl_FID33737.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID326_LC_fresh_R_48_sl_FID33738.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1019_DP\meas_MID327_LC_fresh_R_32_sl_FID33739.dat'

% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID46_LC_R_125_FID30408.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID55_LC_R_75_FID30416.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID67_LC_R_75_60fp_30seg_FID30428.dat'
% measfile = 'C:\Users\Ricardo\Desktop\0928_phantom\meas_MID64_LC_R_125_60fp_FID30425.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID68_LC_offset_R30_without_FID30891.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID70_LC_offset_R60_without_FID30893.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID84_LC_offset_R30_without_125_FID30907.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID88_LC_offset_R60_without_125_FID30911.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID67_LC_offset_R30_with_FID30890.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID69_LC_offset_R60_with_FID30892.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID83_LC_offset_R30_with_125_FID30906.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1001_phan7\meas_MID85_LC_offset_R60_with_125_FID30908.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1004_PHAN8\meas_MID16_LC_fresh_R_FID31585.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1004_PHAN8\meas_MID19_LC_fresh_R_all_yes_FID31588.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID54_LC_fresh_R_offset_64v_FID31662.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID75_LC_fresh_R_offset_64v_SHAX_FID31682.dat'


% measfile = 'C:\Users\Ricardo\Desktop\1012\meas_MID136_LC_fresh_R_64_FID32633.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1012\meas_MID138_LC_fresh_R_64_SL42_FID32635.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1012\meas_MID139_LC_fresh_R_64_SL42_FID32636.dat'

% measfile = 'C:\Users\Ricardo\Desktop\3110_ES\meas_MID213_LC_R_0090180270_48_FID35790.dat'
% measfile ='C:\Users\Ricardo\Desktop\3110_ES\meas_MID214_LC_R_018090270_48_FID35791.dat'
% measfile ='C:\Users\Ricardo\Desktop\3110_ES\meas_MID205_LC_R_018090270_32_FID35782.dat'
% measfile ='C:\Users\Ricardo\Desktop\3110_ES\meas_MID204_LC_R_0090180270_32_FID35781.dat'
% measfile ='C:\Users\Ricardo\Desktop\1012_KR\meas_MID138_LC_fresh_R_64_SL42_FID32635.dat'

% measfile = 'C:\Users\Ricardo\Desktop\meas_MID69_LC_offset_R60_with_FID30892.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID51_LC_R_0090180270_SL1_FID36498.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID52_LC_R_018090270_SL1_FID36499.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID53_LC_R_027090180_SL1_FID36500.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID62_LC_R_0090180270_SL2_FID36509.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1106_JA\meas_MID63_LC_R_018090270_SL2_FID36510.dat'

% measfile = 'C:\Users\Ricardo\Desktop\meas_MID56_3D_radial_cine_LC_FID37504.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1411_MM\meas_MID166_LC_R_0090180270_FID37874.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BALL\meas_MID188_LC_R_GX_COMP_INT_48_16_16_FID39910.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BALL\meas_MID189_LC_R_GX_COMP_INT_64_16_16_FID39911.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BALL\meas_MID190_LC_R_GX_COMP_INT_32_16_16_FID39912.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BALL\meas_MID191_LC_R_GX_COMP_INT_15_15_16_FID39913.dat'


% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID220_LC_R_GX_STD_INT_48_16_16_FID39942.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID221_LC_R_GX_STD_INT_32_16_16_FID39943.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID222_LC_R_GX_STD_INT_64_16_16_FID39944.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID223_LC_R_GX_STD_INT_15_15_16_FID39945.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1130_PHANTOMS\BOTTLES\meas_MID231_LC_R_GX_COMP_INT_15_15_16_FID39953.dat'

% measfile = 'C:\Users\Ricardo\Desktop\1203_MB\meas_MID131_LC_SSFP_48_FID40126.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1203_MB\meas_MID142_LC_SSFP_48_noint_FID40137.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1203_MB\meas_MID144_LC_COMP_SL1_int_48_FID40139.dat'

% measfile ='C:\Users\Ricardo\Desktop\1203_MB\meas_MID138_LC_SSFP_48_sl2_FID40133.dat'
% measfile ='C:\Users\Ricardo\Desktop\1203_MB\meas_MID140_LC_SSFP_48_sl2_noint_FID40135.dat'
% measfile ='C:\Users\Ricardo\Desktop\1203_MB\meas_MID145_LC_COMP_SL2_int_48_FID40140.dat'

% measfile ='C:\Users\Ricardo\Desktop\1205_YR\meas_MID74_LC_SSFP_STD_32_FID40549.dat'
% measfile='C:\Users\Ricardo\Desktop\1205_YR\meas_MID75_LC_SSFP_COMP_32_FID40550.dat'
% measfile ='C:\Users\Ricardo\Desktop\1205_YR\meas_MID78_LC_SSFP_STD_64_FID40553.dat'
% measfile='C:\Users\Ricardo\Desktop\1205_YR\meas_MID80_LC_SSFP_COMP_64_FID40555.dat'
% measfile='C:\Users\Ricardo\Desktop\1205_YR\meas_MID101_LC_SSFP_COMP_48_LV_FID40576.dat'

% measfile='C:\Users\Ricardo\Desktop\0812_bottlephantoms\meas_MID151_LC_R_GX_COMP_INT_15_FID41259.dat'
%meas_MID147_LC_R_GX_COMP_INT_32_FID41255
% meas_MID151_LC_R_GX_COMP_INT_15_FID41259

%READ IN RAW TWIX DATA FILE
dispopt = 'on';
filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);
phases = sMDH.sLC.ushPhase; %PHASES
coils= sMDH.sLC.ulChannelId;%COILS
measure= sMDH.sLC.ushRepetition; %MEASURES / PASSES
part = sMDH.sLC.ushPartition;   %SLICES
Np = sMDH.sLC.ushLine;      %PROJECTIONS PER PASS
Nx = sMDH.ushSamplesInScan; %SAMPLES
%GENERATE PROJECTION ANGLES DEPENDING ON ODD OR EVEN PROJECTIONS
if mod(Np,2) == 0
phival=[0:pi/Np:pi-pi/Np];
end
if mod(Np,2) ~=0
    phival=[0:2*pi/Np:2*pi-2*pi/Np];
end
%TAKE OUTPUT OF READIN FUNCTION AND ORGANISE INTO DESIRED FORMAT
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

%INITIALIZE TO SAVE COMPUTATION TIME
ft_LC_0 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
ft_LC_90 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
ft_LC_180 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
ft_LC_270 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
add = zeros(phases, coils, oversmpl*res,oversmpl*res); 
sub = zeros(phases, coils, oversmpl*res,oversmpl*res); 
four = zeros(phases, coils, oversmpl*res,oversmpl*res); 

%create dcf-matrix (Ram-lak-filter) (same for all coildata)
% dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
% % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
% dcf2=dcf1/max(dcf1);
% dcf=repmat(dcf2,[Np 1]); %384*39 matrix
% figure;plot(dcf(1,:));
% end % if correct_all
% locate kspace samples in complex space (same for all coildata)

%GENERATE TRAJECTORIES FOR GRIDDING
loc=zeros(measure,Nx,Np);
location =zeros(Nx,Np);

for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
for j=[1:Nx] 
    for k=[1:Np] 
         if m == 1
             offset = 0;    %NO OFFSET FOR FIRST PASS
         end
         if m > 1 && mod(Np,2) ==0
             offset = (pi)/(measure*Np);    %SPECIFY OFFSET FOR EVEN Np
         end
         if m >1 && mod(Np,2) ~=0
             offset = (2*pi)/(measure*Np);  %SPECIFY OFFSET FOR ODD Np
         end
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)

        pha=phival(k)+((m-1)*offset);
        loc(m,j,k)=mag*exp(1i*pha);

    end
end
end
%PLOT TRAJECTORIES
% h1 = plot(real(squeeze(loc(1,:,1))),imag(squeeze(loc(1,:,1))), 'color','c');hold on 
% h2 = plot(real(squeeze(loc(2,:,1))),imag(squeeze(loc(2,:,1))), 'color','b');
% h3 = plot(real(squeeze(loc(3,:,1))),imag(squeeze(loc(3,:,1))), 'color','g');
% h4 = plot(real(squeeze(loc(4,:,1))),imag(squeeze(loc(4,:,1))), 'color','m');
% legend
% legend([h1 h2 h3 h4],{'First Pass','Second Pass','Third Pass','Fourth Pass'}, 'Box', 'off')
% h5= plot(real(squeeze(loc(1,:,2:Np))),imag(squeeze(loc(1,:,2:Np))), 'color','c', 'HandleVisibility','off'); 
% h6= plot(real(squeeze(loc(2,:,2:Np))),imag(squeeze(loc(2,:,2:Np))), 'color','b', 'HandleVisibility','off');
% h7= plot(real(squeeze(loc(3,:,2:Np))),imag(squeeze(loc(3,:,2:Np))), 'color','g', 'HandleVisibility','off');
% h8= plot(real(squeeze(loc(4,:,2:Np))),imag(squeeze(loc(4,:,2:Np))), 'color','m', 'HandleVisibility','off');
% set(gca,'Visible','off')
% figure;
% plot(real(squeeze(loc(1,:,:))),imag(squeeze(loc(1,:,:))),'c');hold on
% plot(real(squeeze(loc(2,:,:))),imag(squeeze(loc(2,:,:))),'b')
% plot(real(squeeze(loc(3,:,:))),imag(squeeze(loc(3,:,:))),'g')
% plot(real(squeeze(loc(4,:,:))),imag(squeeze(loc(4,:,:))),'m')
% legend ('First Pass: 0 phase increment','Second Pass: 90 phase increment','Third Pass: 180 phase increment','Fourth Pass: 270 phase increment')


%INITIALIZE
datone=zeros(phases,coils,res,res); %all k-space data 
imone=zeros(phases,coils,res,res);  %all image data 
finalIm=zeros(phases,res,res);   %combined image data from coils 

% FOR EACH SLICE, EACH COIL, EACH PHASE, GRID RADIAL DATA, PERFORM 2D FFT
% AND DEAPODIZATION
for pa=1:part
for cl = 1:coils
    for ph = 1:phases
        size(squeeze(ksampsnewone(pa,cl,ph,:,:)))
        datone(ph,cl,:,:) = gridkb(loc(1,:,:),squeeze(ksampsnewone(pa,cl,ph,:,:)),dcf',oversmpl*res,kwidth,oversmpl);
        imone(ph,cl,:,:) = deapod(fftshift(fft2(fftshift(squeeze(datone(ph,cl,:,:))))),oversmpl,res,kwidth);
        dattwo(ph,cl,:,:) = gridkb(loc(2,:,:),squeeze(ksampsnewtwo(pa,cl,ph,:,:)),dcf',oversmpl*res,kwidth,oversmpl);
        imtwo(ph,cl,:,:) = deapod(fftshift(fft2(fftshift(squeeze(dattwo(ph,cl,:,:))))),oversmpl,res,kwidth);
        datthree(ph,cl,:,:) = gridkb(loc(3,:,:),squeeze(ksampsnewthree(pa,cl,ph,:,:)),dcf',oversmpl*res,kwidth,oversmpl);
        imthree(ph,cl,:,:) = deapod(fftshift(fft2(fftshift(squeeze(datthree(ph,cl,:,:))))),oversmpl,res,kwidth);
        datfour(ph,cl,:,:) = gridkb(loc(4,:,:),squeeze(ksampsnewfour(pa,cl,ph,:,:)),dcf',oversmpl*res,kwidth,oversmpl);
        imfour(ph,cl,:,:) = deapod(fftshift(fft2(fftshift(squeeze(datfour(ph,cl,:,:))))),oversmpl,res,kwidth);

    end
end 
end

%PERFORM LINEAR COMBINATION (CAN ALSO BE DONE IN FOURIER SPACE)
%NOTE WHICH IS 90 AND WHICH IS 180!! DEPENDS ON ACQUISITION ORDER
%FOR STANDARD ORDER FT_LC_90 = IMTWO, FT_LC_180 = IMTHREE
%FOR COMPLEMENTARY ORDER FT_LC_90 = IMTHREE, FT_LC_180 = IMTWO
for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_90(phase, coil,:,:) = (imtwo(phase, coil, :, :)); %CHANGE TO 180 FOR COMPLEMENTARY ORDER
        ft_LC_180(phase, coil,:,:) = (imthree(phase, coil, :, :)); %CHANGE TO 90 FOR COMPLEMENTARY ORDER
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
    end
end

%CREATE AVI FILE
% vidObj = VideoWriter('TRAIL2');
% vidObj.FrameRate = 2;
% open(vidObj);

%PLOT FOUR PASS COMBINATION, AND EACH PHASE CYCLING IMAGE FOR ALL PHASES
oot=1; %COUNTER
while oot<6
for phase = 1:phases

% % 
if oot == 1;
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,104:519,104:519))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = flip(final_img(phase,:,:),3);
    colormap(gray); imagesc(squeeze(final(phase,:,:))); axis equal; caxis([0,max(max(max(final)))]); title('FOUR COMP INT 15p');
%     set(gca,'Xdir','reverse')
end
% % % %  
% % % % %     
if oot==2;
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(phase,:,104:519,104:519))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = (flip(final_img(phase,:,:),3));
colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('180 COMP INT 15p');  
% set(gca,'Xdir','reverse')
end
if oot==3;
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_0(phase,:,104:519,104:519))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = flip(final_img(phase,:,:),3);
colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('0 COMP INT 15p');  
% set(gca,'Xdir','reverse')
end
if oot==4;
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_90(phase,:,104:519,104:519))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = flip(final_img(phase,:,:),3);
colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('90 COMP INT 15p');  
% set(gca,'Xdir','reverse')
end
if oot==5;
figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_270(phase,:,104:519,104:519))));
final_img(phase,:,:) = ((((finalIm(phase,:,:)))));
final(phase,:,:) = flip(final_img(phase,:,:),3);
colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('270 COMP INT 15p');  
% set(gca,'Xdir','reverse')
end
    
%         currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);

%GENERATE NIFTI FILE
voxel = [2 2 10];
myimage=(permute(final,[3 2 1])); % create myimage\
% myimage = flip(myimage, 2);

if oot == 1;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms\','LC_SSFP_COMP_INT_15_FOUR'); end %create a %%file name
if oot == 2;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms\','LC_SSFP_COMP_INT_15_180'); end
if oot ==3;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms\','LC_SSFP_COMP_INT_15_0'); end
if oot == 4;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms\','LC_SSFP_COMP_INT_15_90'); end
if oot == 5;imfilename=sprintf('%s%s','C:\Users\Ricardo\Desktop\0812_bottlephantoms\','LC_SSFP_COMP_INT_15_270'); end
nii= make_nii(myimage, voxel);         
save_nii(nii, imfilename);

%   view_nii(nii);
end
  oot = oot+1;
    clear nii
  end
% close(vidObj);
% 
% 



%%PLOT MORE THAN ONE VOLUNTEER ON SAME FIGURE
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