% 2022/3/6, jie xiang @yale mrrc
% read GRE field mapping

% close all
% clear
% read rawdata
% measfile1 = 'F:\CVCoding\gridding_lcssfp\MeasData\20220106_Phantom_HomoBall\meas_MID00070_FID136157_AdjGre.dat';
measfile1 = 'F:\CVCoding\gridding_lcssfp\MeasData\20220303_SupumLee\Cartesian\meas_MID00150_FID140709_AdjGre.dat';
dispopt = 'on';
[ MultiRAID1, measfile1, ~ ] = ReadSiemensMeasVD13_idea_memoryMod(measfile1, dispopt);
rawdata0 = MultiRAID1(1).rawdata;

% rawdata1 = zeros(2,88,128);
% rawdata2 = zeros(2,88,128);
% recon the first imag, channel 1
rawdata1(1,:,:) = rawdata0(3:4:end,:);
imag1(1,:,:)=fft2c(fftshift(squeeze(rawdata1(1,:,:))));
% channel 2
rawdata1(2,:,:) = rawdata0(4:4:end,:);
imag1(2,:,:)=fft2c(fftshift(squeeze(rawdata1(1,:,:))));

img1 = AdaptiveCombine(imag1);
figure,imshow(abs(img1),[]);

% recon the second imag
rawdata2(1,:,:) = rawdata0(5:4:end,:);
imag2(1,:,:)=fft2c(fftshift(squeeze(rawdata2(1,:,:))));
% channel 2
rawdata2(2,:,:) = rawdata0(6:4:end,:);
imag2(2,:,:)=fft2c(fftshift(squeeze(rawdata2(1,:,:))));

img2 = AdaptiveCombine(imag2);
figure,imshow(abs(img2),[]);

% generate omega map
TE1 = 2.98; TE2 = 4.91; Del_TE = TE2 - TE1;
OmegaMap = angle(img1./img2)./Del_TE;
figure, imshow(OmegaMap,[])

roi = images.roi.AssistedFreehand;
draw(roi);
mask = createMask(roi);
OmegaMap3 = OmegaMap.*mask;
figure, imshow(OmegaMap3.*1000,[]),colorbar,title('Omega Map,Hz')
