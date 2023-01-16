%RECONSTRUCT WITH WEIGHTED COMBINATION METHOD FOR LC-SSFP
%IDENTICAL METHOD UNTIL COMBINATION STEP

%set initial parameters
res=160;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\0914_YR2twix');
folder = 'C:\Users\Ricardo\Desktop\0914_YR2twix';
filepath = '';

%SPECIFY INPUT FILE
% measfile = 'C:\Users\Ricardo\Desktop\1025_BR\timtriob_20181025_064338988\meas_MID256_LC_R_0090180270_FID34977.dat'

%READ IN TWIX FILE AND ORGANISE
dispopt = 'on';
filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);
phases = sMDH.sLC.ushPhase;
coils=sMDH.sLC.ulChannelId; %sMDH.ushUsedChannels;
measure= sMDH.sLC.ushRepetition;
Np = sMDH.sLC.ushLine;
Nx = sMDH.ushSamplesInScan;
if mod(Np,2) == 0
phival=[0:pi/Np:pi-pi/Np];
end
if mod(Np,2) ~=0
    phival=[0:2*pi/Np:2*pi-2*pi/Np];
end
nunder=1;   % undersampling factor
phival= phival(1:nunder:end);
ksampsone = first;
ksampsnewone=permute(ksampsone,[1 2 4 3]);
ksampsnewone= ksampsnewone(:,:,1:nunder:end,:);
ksampstwo = second;
ksampsnewtwo=permute(ksampstwo,[1 2 4 3]);
ksampsnewtwo= ksampsnewtwo(:,:,1:nunder:end,:);
ksampsthree = third;
ksampsnewthree=permute(ksampsthree,[1 2 4 3]);
ksampsnewthree= ksampsnewthree(:,:,1:nunder:end,:);
ksampsfour = fourth;
ksampsnewfour=permute(ksampsfour,[1 2 4 3]);
ksampsnewfour= ksampsnewfour(:,:,1:nunder:end,:);

interleaving = 1; %or 2 for off
order = 1; % standard = 1, complementary = 2

% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
k=[-.5:1/Nx:.5-1/Nx];
dcf1=abs(k).*sinc(k);
dcf=repmat(dcf1,[Np 1]); 

%INITIALIZE
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

%GENERATE RADIAL TRAJECTIORIES FOR GRIDDING
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

% figure;
%perform gridding (different for all coildata) and create image with
%fourier transform
% datone=zeros(phases,coils,res,res); %all k-space data (empty matrix)
% imone=zeros(phases,coils,res,res);  %all image data (empty matrix)
% finalIm=zeros(phases,res,res);   %combined image data from coils (empty matrix)

%PERFORM GRIDDING AND 2D FFT
for ph=1:phases
    for cl=1:coils
        size(squeeze(ksampsnewone(cl,ph,:,:)))
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

%THIS STEP SLIGHTLY DIFFERS FROM NORMAL LC-SSFP
% P IS THE WEIGHTING FACTOR
% EACH COMBINED IMAGE IS PLOTTED USING THE SPECIFIED P VALUES 
% REFERENCE: CUKUR ET AL. ENHANCED SPECTRAL SHAPING IN SSFP IMAGING
for p = [-0.5 0 0.5 1 2 3 4 5 10]
%p = -1;
for coil = 1:coils
    for phase = 1:phases
        if order == 1
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_90(phase, coil,:,:) = (imtwo(phase, coil, :, :));
        ft_LC_180(phase, coil,:,:) = (imthree(phase, coil, :, :));
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        add(phase, coil,:,:) = (abs((abs(ft_LC_0(phase, coil,:,:)).^p).*(ft_LC_0(phase, coil,:,:))+((sqrt(-1))*((abs(ft_LC_270(phase, coil, :,:)).^p).*(ft_LC_270(phase, coil, :,:))))).^(1/(1+p)));
        sub(phase, coil,:,:) = (abs((abs(ft_LC_0(phase, coil,:,:)).^p).*(ft_LC_0(phase, coil,:,:))-((sqrt(-1))*((abs(ft_LC_180(phase, coil, :,:)).^p).*(ft_LC_180(phase, coil, :,:))))).^(1/(1+p)));
        four(phase, coil,:,:) = (abs(((abs(ft_LC_0(phase, coil,:,:))).^p.*(ft_LC_0(phase, coil,:,:)))+((abs(ft_LC_90(phase, coil,:,:))).^p.*(ft_LC_90(phase, coil,:,:)))+((abs(ft_LC_180(phase, coil,:,:))).^p.*(ft_LC_180(phase, coil,:,:)))+((abs(ft_LC_270(phase, coil,:,:))).^p.*(ft_LC_270(phase, coil,:,:))))).^(1/(1+p)); 
        end
        if order == 2 
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_180(phase, coil,:,:) = (imtwo(phase, coil, :, :));
        ft_LC_90(phase, coil,:,:) = (imthree(phase, coil, :, :));
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        add(phase, coil,:,:) = (abs((abs(ft_LC_0(phase, coil,:,:)).^p).*(ft_LC_0(phase, coil,:,:))+((sqrt(-1))*((abs(ft_LC_270(phase, coil, :,:)).^p).*(ft_LC_270(phase, coil, :,:))))).^(1/(1+p)));
        sub(phase, coil,:,:) = (abs((abs(ft_LC_0(phase, coil,:,:)).^p).*(ft_LC_0(phase, coil,:,:))-((sqrt(-1))*((abs(ft_LC_180(phase, coil, :,:)).^p).*(ft_LC_180(phase, coil, :,:))))).^(1/(1+p)));
        four(phase, coil,:,:) = (abs(((abs(ft_LC_0(phase, coil,:,:))).^p.*(ft_LC_0(phase, coil,:,:)))+((abs(ft_LC_90(phase, coil,:,:))).^p.*(ft_LC_90(phase, coil,:,:)))+((abs(ft_LC_180(phase, coil,:,:))).^p.*(ft_LC_180(phase, coil,:,:)))+((abs(ft_LC_270(phase, coil,:,:))).^p.*(ft_LC_270(phase, coil,:,:))))).^(1/(1+p));
        end
end

%CREATE AVI FILE
% vidObj = VideoWriter(char(sprintf( "%s%s%d", 'DP_090180270  ', 'p = ', p)));
% open(vidObj);
for phase = 1:phases
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(add(phase,:,:,:))));
%     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD Phase   ');
% % % 
%  figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(add(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'ADD  ', 'p = ', p));
% % % % 
% % % % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(sub(phase,:,:,:))));
% % % % % % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('SUBTRACT Phase  ');
% % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(sub(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'SUBTRACT  ', 'p = ', p));
% % % % % % %     
% % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(four(phase,:,:,:))));
% % %     colormap(gray); imagesc(fliplr(((squeeze(finalIm(phase,:,:))))); axis equal; title('FOUR Phase ) ');

figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
    colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'FOUR  ', 'p = ', p));
 
% %     
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title('180 -OFFSET');  
    
%         currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
end

% close(vidObj);
end