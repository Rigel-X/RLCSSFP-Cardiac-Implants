% gridding based on hargreaves spiral gridding code.

%set initial parameters
res=160;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\0914_YR2twix');
folder = 'C:\Users\Ricardo\Desktop\0914_YR2twix';
filepath = '';

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

% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\meas_MID239_LC_fresh_R_32_sl_FID33651.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\meas_MID241_LC_fresh_R_48_sl_FID33653.dat'
% measfile = 'C:\Users\Ricardo\Desktop\1018_IO\meas_MID243_LC_fresh_R_64_sl_FID33655.dat'

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


dispopt = 'on';
filename = 'meas_MID138_LC_R_180_FID29431.dat';   % radial 128 spoke, long TR, 100 mm FOV
[first, second, third, fourth, rawdata, newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17measures(measfile, dispopt);
phases = sMDH.sLC.ushPhase;
part = sMDH.sLC.ushPartition;
coils= sMDH.sLC.ulChannelId; %sMDH.sLC.ushPartition; sMDH.ushUsedChannels;
measure= sMDH.sLC.ushRepetition;
Np = sMDH.sLC.ushLine;
Nx = sMDH.ushSamplesInScan;
phival=[0:pi/Np:pi-pi/Np];
nunder=1;% undersampling factor
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

% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
k=[-.5:1/Nx:.5-1/Nx];
dcf1=abs(k).*sinc(k);
dcf=repmat(dcf1,[Np 1]); 

ft_LC_0 = zeros(part, coils, res, res); 
ft_LC_90 = zeros(part, coils, res, res); 
ft_LC_180 = zeros(part, coils, res, res); 
ft_LC_270 = zeros(part, coils, res, res); 
add = zeros(part, coils, res, res); 
sub = zeros(part, coils, res, res); 
four = zeros(part, coils, res, res); 

%create dcf-matrix (Ram-lak-filter) (same for all coildata)
% dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
% % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
% dcf2=dcf1/max(dcf1);
% dcf=repmat(dcf2,[Np 1]); %384*39 matrix
% figure;plot(dcf(1,:));

%end % if correct_all
%locate kspace samples in complex space (same for all coildata)
loc=zeros(measure,Nx,Np);

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
        loc(m,j,k)=mag*exp(i*pha);

    end
end
% plot(real(squeeze(loc(m,:,:))),imag(squeeze(loc(m,:,:))), 'c')
% hold on 
end

% figure;
%perform gridding (different for all coildata) and create image with
%fourier transform
% datone=zeros(phases,coils,res,res); %all k-space data (empty matrix)
% imone=zeros(phases,coils,res,res);  %all image data (empty matrix)
% finalIm=zeros(phases,res,res);   %combined image data from coils (empty matrix)

for pa=1:part
    for cl=1:coils
        size(squeeze(ksampsnewone(cl,pa,:,:)))
        datone(pa,cl,:,:) = gridkb(loc(1,:,:),squeeze(ksampsnewone(cl,pa,:,:)),dcf',res,kwidth,oversmpl);
        imone(pa,cl,:,:) = (fft2(squeeze(datone(pa,cl,:,:))));
        dattwo(pa,cl,:,:) = gridkb(loc(2,:,:),squeeze(ksampsnewtwo(cl,pa,:,:)),dcf',res,kwidth,oversmpl);
        imtwo(pa,cl,:,:) = (fft2(squeeze(dattwo(pa,cl,:,:))));
        datthree(pa,cl,:,:) = gridkb(loc(3,:,:),squeeze(ksampsnewthree(cl,pa,:,:)),dcf',res,kwidth,oversmpl);
        imthree(pa,cl,:,:) = (fft2(squeeze(datthree(pa,cl,:,:))));
        datfour(pa,cl,:,:) = gridkb(loc(4,:,:),squeeze(ksampsnewfour(cl,pa,:,:)),dcf',res,kwidth,oversmpl);
        imfour(pa,cl,:,:) = (fft2(squeeze(datfour(pa,cl,:,:))));
    end
%     finalIm(ph,:,:)=mergeCoilImages(squeeze(abs(im(ph,:,:,:))));
%     colormap(gray); imagesc(fftshift(squeeze(finalIm(ph,:,:))));
   
end 

for coil = 1:coils
    for pa = 1:part
        ft_LC_0(pa, coil,:,:) = (imone(pa, coil, :, :));
        ft_LC_270(pa, coil,:,:) = (imtwo(pa, coil, :, :));
        ft_LC_90(pa, coil,:,:) = (imthree(pa, coil, :, :));
        ft_LC_180(pa, coil,:,:) = (imfour(pa, coil, :, :));
        add(pa, coil,:,:) = ft_LC_0(pa, coil,:,:)+((sqrt(-1))*(ft_LC_180(pa, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_90(phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
        sub(pa, coil,:,:) = ft_LC_0(pa, coil,:,:)-((sqrt(-1))*(ft_LC_180(pa, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_90(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
        four(pa, coil,:,:) = (ft_LC_0(pa, coil,:,:)) + (ft_LC_180(pa, coil, :,:)) + (ft_LC_270(pa, coil,:,:)) + (ft_LC_90(pa, coil, :,:));
    end
end

% vidObj = VideoWriter('BR_027090180_270');
% open(vidObj);
% phases = 1;
% phases = 1;
% for coil = 1:coils
for pa = 1:part
% % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(add(phase,:,:,:))));
% %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('ADD Phase   ');
% % % % % % % % 
%  figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(add(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title('ADD   ');
% % % % % 
% % % % % % % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(sub(phase,:,:,:))));
% % % % % % % % % %     colormap(gray); imagesc((squeeze(finalIm(phase,:,:)))); axis equal; title('SUBTRACT Phase  ');
% % % % % % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(sub(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title('SUBTRACT   ');
% % % % % % % % %     
% % % % % % figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(four(phase,:,:,:))));
% % % % % %     colormap(gray); imagesc(fliplr(((squeeze(finalIm(phase,:,:))))); axis equal; title('FOUR Phase ) ');
% % 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(four(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(phase,:,:))))); axis equal; title('FOUR   ');
 
% % %     
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(abs(ft_LC_180(pa,:,:,:))));
%     colormap(gray); imagesc(fliplr(fftshift(squeeze(finalIm(pa,:,:))))); axis equal; title('180 ');  
    
    figure; colormap(gray);imagesc(fftshift(squeeze(abs((ft_LC_180(pa,coil,:,:)))))); axis equal;title('180');
%     
%         currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
end
% end
% close(vidObj);