clear all

%SET INITIAL PARAMTERS 
res=192;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\');
folder = 'C:\Users\Ricardo\Desktop\';
filepath = '';
close all;

%%
%%% CHenxi HU>>>>>>>>>>>>>>>>>
%%%  CHANGE THE DIRECTORY/ADD iN STANDARD ORDER FILES >>>>>>>>>>>>>>>>>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 measfile0{1} = 'C:\Users\Ricardo\Desktop\JaimesFiles\Phantoms\0812_bottlephantoms\meas_MID136_LC_R_GX_COMP_INT_48_FID41244.dat' %Chenxi
 measfile0{2} = 'C:\Users\Ricardo\Desktop\JaimesFiles\Phantoms\0812_bottlephantoms\meas_MID138_LC_R_GX_COMP_NOINT_48_FID41246.dat' %chenxi
%%%%%%%%%%%%%%%%%%%%%%%%%%
% % measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID256_LC_R_0090180270_FID34977.dat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 20MM SLICE %%%

%measfile0{1} = 'C:\Users\Ricardo\Desktop\JaimesFiles\phantom02132019\20mmST\meas_MID243_LC_SSFP_COMP20MM_16_60_FID54411.dat' %Chenxi
%measfile0{2} = 'C:\Users\Ricardo\Desktop\JaimesFiles\phantom02132019\20mmST\meas_MID254_LC_SSFP_COMP20MM_NONINT_16_60_FID54422.dat' %chenxi


%%% 4MM SLICE %%%
%measfile0{1} = 'C:\Users\Ricardo\Desktop\JaimesFiles\phantom02132019\thin\meas_MID219_LC_SSFP_COMP_384_60_FID54387.dat' %Chenxi
%measfile0{2} = 'C:\Users\Ricardo\Desktop\JaimesFiles\phantom02132019\thin\meas_MID230_LC_SSFP_COMP_NONINT_384_60_FID54398.dat' %chenxi
%%
combMethod = 'CS-SSFP';
% combMethod = 'SOS-SSFP';
% combMethod = 'MS-SSFP';
for measfileIdx = 2
    measfile = measfile0{measfileIdx};

%% READ IN RAW TWIX DATA FILE
dispopt = 'on';
% filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17_2DLCSSFP(measfile, dispopt);
%%
for Np = [32,48,64,96,128,192,384]
    interval = 384/Np;
    first1 = first(:,:,1:interval:end,:);
    second1 = second(:,:,1:interval:end,:);
    third1 = third(:,:,1:interval:end,:);
    fourth1 = fourth(:,:,1:interval:end,:);

% crop k-space
Nx = sMDH.ushSamplesInScan; %SAMPLES
% cropind = Nx/2-Nx/8:Nx/2+Nx/8-1;
% Nx = length(cropind);
% first = first(:,:,:,cropind);
% second = second(:,:,:,cropind);
% third = third(:,:,:,cropind);
% fourth = fourth(:,:,:,cropind);

phases = sMDH.sLC.ushPhase; %PHASES
coils= sMDH.sLC.ulChannelId;%COILS
measure= sMDH.sLC.ushRepetition; %MEASURES / PASSES
part = sMDH.sLC.ushPartition;   %SLICES
% Np = sMDH.sLC.ushLine;      %PROJECTIONS PER PASS

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
ksampsone = first1;
ksampsnewone=permute(ksampsone,[1 2 4 3]);
ksampsnewone= ksampsnewone(:,:,:,:,1:nunder:end);
ksampstwo = second1;
ksampsnewtwo=permute(ksampstwo,[1 2 4 3]);
ksampsnewtwo= ksampsnewtwo(:,:,:,:,1:nunder:end);
ksampsthree = third1;
ksampsnewthree=permute(ksampsthree,[1 2 4 3]);
ksampsnewthree= ksampsnewthree(:,:,:,:,1:nunder:end);
ksampsfour = fourth1;
ksampsnewfour=permute(ksampsfour,[1 2 4 3]);
ksampsnewfour= ksampsnewfour(:,:,:,:,1:nunder:end);

%IS DATA INTERLEAVED OR NOT?
if mod(measfileIdx,2)==1
    interleaving = 1; %1 = on, 2 = off
else
    interleaving = 2;
end
%WHICH ORDER OF ACQUISITION WAS USED
order = 2; % standard order = 1, complementary order = 2
res=Nx/2;
%% reconstruction
% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
shift=0.0;
k=[-.5-shift:1/Nx:.5-1/Nx-shift];
dcf1=abs(k);%.*cos(pi*k);
dcf=repmat(dcf1,[Np 1]); 

%INITIALIZE TO SAVE COMPUTATION TIME
ft_LC_0 = zeros(phases, coils, res,res); 
ft_LC_90 = zeros(phases, coils, res,res); 
ft_LC_180 = zeros(phases, coils, res,res); 
ft_LC_270 = zeros(phases, coils, res,res); 
add = zeros(phases, coils, res,res); 
sub = zeros(phases, coils, res,res); 
four = zeros(phases, coils, res,res); 

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

if interleaving == 1
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
end

if interleaving==2;
    for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
for j=[1:Nx] 
    for k=[1:Np] 

        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)

        pha=phival(k);
        loc(m,j,k)=mag*exp(1i*pha);

    end
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
% datone=zeros(phases,coils,res,res); %all k-space data 
% imone=zeros(phases,coils,res,res);  %all image data 
% finalIm=zeros(phases,res,res);   %combined image data from coils 

% FOR EACH SLICE, EACH COIL, EACH PHASE, GRID RADIAL DATA, PERFORM 2D FFT
% AND DEAPODIZATION
% for pa=1:part
%%
clear datone imone dattwo imtwo datthree imthree datfour imfour
warning('off');
for cl = 1:coils
    for ph = 1:phases
%         size(squeeze(ksampsnewone(cl,ph,:,:)))
        datone(ph,cl,:,:) = gridkb(loc(1,:,:),squeeze((ksampsnewone(cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        imone(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datone(ph,cl,:,:)))))));
%         imone(ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(datone(ph,cl,:,:))))),oversmpl,res,kwidth);
        dattwo(ph,cl,:,:) = gridkb(loc(2,:,:),squeeze((ksampsnewtwo(cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        imtwo(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(dattwo(ph,cl,:,:)))))));
        datthree(ph,cl,:,:) = gridkb(loc(3,:,:),squeeze((ksampsnewthree(cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        imthree(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datthree(ph,cl,:,:)))))));
        datfour(ph,cl,:,:) = gridkb(loc(4,:,:),squeeze((ksampsnewfour(cl,ph,:,:))),dcf',res,kwidth,oversmpl);
        imfour(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datfour(ph,cl,:,:)))))));

    end
end 
warning('on');
% end
%%
%PERFORM LINEAR COMBINATION (CAN ALSO BE DONE IN FOURIER SPACE)
%NOTE WHICH IS 90 AND WHICH IS 180!! DEPENDS ON ACQUISITION ORDER
%FOR STANDARD ORDER FT_LC_90 = IMTWO, FT_LC_180 = IMTHREE
%FOR COMPLEMENTARY ORDER FT_LC_90 = IMTHREE, FT_LC_180 = IMTWO
if order == 1
for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_90(phase, coil,:,:) = (imtwo(phase, coil, :, :)); 
        ft_LC_180(phase, coil,:,:) = (imthree(phase, coil, :, :)); 
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        % deapodization
        temp = squeeze(ft_LC_0(phase, coil,:,:));
        ft_LC_0(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_90(phase, coil,:,:));
        ft_LC_90(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_180(phase, coil,:,:));
        ft_LC_180(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_270(phase, coil,:,:));
        ft_LC_270(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        
        
        add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        if strcmp(combMethod, 'CS-SSFP')
            four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
        elseif strcmp(combMethod, 'SOS-SSFP') %SOS-SSFP
            four(phase, coil,:,:) = sqrt(abs(ft_LC_0(phase, coil,:,:)).^2 + abs(ft_LC_90(phase, coil, :,:)).^2 + abs(ft_LC_180( phase, coil,:,:)).^2 + abs(ft_LC_270( phase, coil, :,:)).^2);
        else
            four(phase, coil,:,:) = (abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:)));
        end
    end
end
end
if order == 2
    for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_180(phase, coil,:,:) = (imtwo(phase, coil, :, :)); 
        ft_LC_90(phase, coil,:,:) = (imthree(phase, coil, :, :));
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        
        % deapodization
        temp = squeeze(ft_LC_0(phase, coil,:,:));
        ft_LC_0(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_90(phase, coil,:,:));
        ft_LC_90(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_180(phase, coil,:,:));
        ft_LC_180(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        temp = squeeze(ft_LC_270(phase, coil,:,:));
        ft_LC_270(phase, coil,:,:) = deapod( temp,oversmpl,res, kwidth);
        
%         if coil<2
%         figure(1); hold on;
%         imagesc(angle(squeeze(ft_LC_0(1,1,:,:))));title('1'); colormap(jet);pause(0.5);
%         imagesc(angle(squeeze(ft_LC_90(1,1,:,:))));title('2'); colormap(jet);pause(0.5);
%         imagesc(angle(squeeze(ft_LC_180(1,1,:,:))));title('3'); colormap(jet);pause(0.5);
%         imagesc(angle(squeeze(ft_LC_270(1,1,:,:))));title('4'); colormap(jet);pause(0.5);
%         end
%         
        add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        if strcmp(combMethod, 'CS-SSFP')
            four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
        elseif strcmp(combMethod, 'SOS-SSFP') %SOS-SSFP
            four(phase, coil,:,:) = sqrt(abs(ft_LC_0(phase, coil,:,:)).^2 + abs(ft_LC_90(phase, coil, :,:)).^2 + abs(ft_LC_180( phase, coil,:,:)).^2 + abs(ft_LC_270( phase, coil, :,:)).^2);
        else
            four(phase, coil,:,:) = (abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:)));
        end
    end
end
end
%
four2 = sqrt(sum(abs(four).^2,2));
four2 = squeeze(four2);
four2 = permute(four2,[2,3,1]);
% close all;
figure; hold on;
for idx = 1:phases
    imagesc(abs(four2(:,:,idx)));colormap(gray);pause(0.1)
end
%%
%CREATE AVI FILE
% vidObj = VideoWriter('TRAIL2');
% vidObj.FrameRate = 2;
% open(vidObj);
clear finalIm final_img final
siize = res;
start = 1;
finito = res;
start = 1;
finito = siize;
%PLOT FOUR PASS COMBINATION, AND EACH PHASE CYCLING IMAGE FOR ALL PHASES
oot=1; %COUNTER
while oot<6
    for phase = 1:3%phases

        % % 
        if oot == 1;
        %figure;
        finalIm(phase,1,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,start:finito, start:finito))));
        final_img(phase,1,:,:) = ((((finalIm(phase,1,:,:)))));
        final(phase,1,:,:) = flip(final_img(phase,1,:,:),4);
        
        finalIm(phase,2,:,:)=2*pi*double(squeeze(angle(four(phase,1,start:finito, start:finito)))<0)+(squeeze(angle(four(phase,1,start:finito, start:finito))));
        final_img(phase,2,:,:) = ((((finalIm(phase,2,:,:)))));
        final(phase,2,:,:) = flip(final_img(phase,2,:,:),4);
          %  colormap(gray); imagesc(squeeze(final(phase,:,:))); axis equal; caxis([0,max(max(max(final)))]); title('FOUR STD INT 48p');
        %     set(gca,'Xdir','reverse')
        end
        % % % %  
        % % % % %     
        if oot==2;
        %figure;
        finalIm(phase,1,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(phase,:,start:finito, start:finito))));
        final_img(phase,1,:,:) = ((((finalIm(phase,1,:,:)))));
        final(phase,1,:,:) = (flip(final_img(phase,1,:,:),4));
        
        finalIm(phase,2,:,:)=2*pi*double(squeeze(angle(ft_LC_180(phase,1,start:finito, start:finito)))<0)+(squeeze(angle(ft_LC_180(phase,1,start:finito, start:finito))));
        final_img(phase,2,:,:) = ((((finalIm(phase,2,:,:)))));
        final(phase,2,:,:) = (flip(final_img(phase,2,:,:),4));
        %colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('180 STD INT 48p');  
        % set(gca,'Xdir','reverse')
        end
        if oot==3;
        %figure;
        finalIm(phase,1,:,:)=mergeCoilImages(squeeze(abs(ft_LC_0(phase,:,start:finito, start:finito))));
        final_img(phase,1,:,:) = ((((finalIm(phase,1,:,:)))));
        final(phase,1,:,:) = flip(final_img(phase,1,:,:),4);
        
        finalIm(phase,2,:,:)=2*pi*double(squeeze(angle(ft_LC_0(phase,1,start:finito, start:finito)))<0)+(squeeze(angle(ft_LC_0(phase,1,start:finito, start:finito))));
        final_img(phase,2,:,:) = ((((finalIm(phase,2,:,:)))));
        final(phase,2,:,:) = flip(final_img(phase,2,:,:),4);
        %colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('0 STD INT 48p');  
        % set(gca,'Xdir','reverse')
        end
        if oot==4;
        %figure;
        finalIm(phase,1,:,:)=mergeCoilImages(squeeze(abs(ft_LC_90(phase,:,start:finito, start:finito))));
        final_img(phase,1,:,:) = ((((finalIm(phase,1,:,:)))));
        final(phase,1,:,:) = flip(final_img(phase,1,:,:),4);
        
        finalIm(phase,2,:,:)=2*pi*double(squeeze(angle(ft_LC_90(phase,1,start:finito, start:finito)))<0)+(squeeze(angle(ft_LC_90(phase,1,start:finito, start:finito))));
        final_img(phase,2,:,:) = ((((finalIm(phase,2,:,:)))));
        final(phase,2,:,:) = flip(final_img(phase,2,:,:),4);
        %colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('90 STD INT 48p');  
        % set(gca,'Xdir','reverse')
        end
        if oot==5;
        %figure;
        finalIm(phase,1,:,:)=mergeCoilImages(squeeze(abs(ft_LC_270(phase,:,start:finito, start:finito))));
        final_img(phase,1,:,:) = ((((finalIm(phase,1,:,:)))));
        final(phase,1,:,:) = flip(final_img(phase,1,:,:),4);
        
        finalIm(phase,2,:,:)=2*pi*double(squeeze(angle(ft_LC_270(phase,1,start:finito, start:finito)))<0)+(squeeze(angle(ft_LC_270(phase,1,start:finito, start:finito))));
        final_img(phase,2,:,:) = ((((finalIm(phase,2,:,:)))));
        final(phase,2,:,:) = flip(final_img(phase,2,:,:),4);
        %colormap(gray); imagesc((squeeze((final(phase,:,:))))); axis equal; caxis([0,max(max(max(final)))]); title('270 STD INT 48p');  
        % set(gca,'Xdir','reverse')
        end

        %         currFrame = getframe(gcf);
        %     writeVideo(vidObj,currFrame);

        %GENERATE NIFTI FILE
        voxel = [2 2 10];
        myimage=(permute(final,[4 3 1 2])); % create myimage\
        % myimage = flip(myimage, 2);
        if oot == 1, imgfour = myimage(:,:,1,1); end
        if oot == 2, imgssfp = myimage(:,:,1,1); end

        if oot == 1;imfilename=sprintf('%s%s%s%s','C:\Users\Ricardo\Desktop\temp\',measfile(42:52), combMethod, [measfile(end-21:end-4) '_FOUR.nii']); end %create a %%file name
        if oot == 2;imfilename=sprintf('%s%s%s%s','C:\Users\Ricardo\Desktop\temp\',measfile(42:52),combMethod,[measfile(end-21:end-4) '_180.nii']); end
        if oot ==3;imfilename=sprintf('%s%s%s%s','C:\Users\Ricardo\Desktop\temp\',measfile(42:52),combMethod,[measfile(end-21:end-4) '_0.nii']); end
        if oot == 4;imfilename=sprintf('%s%s%s%s','C:\Users\Ricardo\Desktop\temp\',measfile(42:52),combMethod,[measfile(end-21:end-4) '_90.nii']); end
        if oot == 5;imfilename=sprintf('%s%s%s%s','C:\Users\Ricardo\Desktop\temp\',measfile(42:52),combMethod,[measfile(end-21:end-4) '_270.nii']); end
        nii= make_nii(myimage(:,:,:,1), voxel);         
        if oot==1, save_nii(nii, imfilename);  end

    %   view_nii(nii);
    end
    oot = oot+1;
    clear nii
end

save(['C:\Users\Ricardo\Desktop\temp\matfile\' measfile(42:52) combMethod measfile(end-21:end-4) '_Np' num2str(Np) '.mat'], 'imgfour', 'imgssfp', 'four', 'four2');
% close(vidObj);
% 
% 
end
end