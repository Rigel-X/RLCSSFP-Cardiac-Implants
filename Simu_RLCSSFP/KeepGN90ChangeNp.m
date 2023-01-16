% simulated LC-SSFP with Gaussian Noise 90 and different
% projections Np, all recon data under folder ComplexSumGN90ChangeNp\
% plot figures used in JMRI revised manuscript

% 2022/5/31, jie xiang @yale mrrc


close all
namelist = dir('F:\CVCoding\gridding_lcssfp\Simulation\ComplexSumGN90ChangeNp\*.mat');
addpath 'F:\CVCoding\gridding_lcssfp'

for NoProject = 1:length(namelist)
filename = namelist(NoProject).name;
load(filename)
% cd('F:\CVCoding\gridding_lcssfp')
RAW1 = Ksp1;
RAW2 = Ksp2;
RAW3 = Ksp3;
RAW4 = Ksp4;
%% Image with no noise
[ndy,Nx,Nspoke] = size(Ksp1);
isaddnoise = 0;
if isaddnoise == 1
    GNnoiseSNR = 100;
    for np = 1:20
    for i=1:4
        for j = 1:Nspoke
            Ksp1(i,:,j,np) = awgn(squeeze(RAW1(i,:,j)),GNnoiseSNR);
            Ksp2(i,:,j,np) = awgn(squeeze(RAW2(i,:,j)),GNnoiseSNR);
            Ksp3(i,:,j,np) = awgn(squeeze(RAW3(i,:,j)),GNnoiseSNR);
            Ksp4(i,:,j,np) = awgn(squeeze(RAW4(i,:,j)),GNnoiseSNR);
        end
    end
    end
else
    for np = 1:20
        Ksp1(:,:,:,np) = RAW1;
        Ksp2(:,:,:,np) = RAW2;
        Ksp3(:,:,:,np) = RAW3;
        Ksp4(:,:,:,np) = RAW4;
    end
end
% figure,plot(abs(squeeze(Ksp1(1,:,1))))
% figure,plot(abs(squeeze(RAW1(1,:,1))))
% figure,plot(abs(squeeze(RAW1(1,:,1)))-abs(squeeze(Ksp1(1,:,1))))

%% Streak Artifact Ratio
Nx = size(lcSSFP,1);
FilterBlock = zeros(Nx,Nx);
FilterBlock((Nx/2-30):(Nx/2+30),(Nx/2-30):(Nx/2+30)) = 1;
for loop = 1:3
for np = 1:20
if loop == 1
    Ij = squeeze(abs(lcSSFP(:,:,np)));
elseif loop == 2
    Ij = squeeze(abs(idSSFP(:,:,np)));
else
    Ij = squeeze(abs(ipSSFP(:,:,np)));
end
IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
Lj = fftshift(fft2(ifftshift(IjTemp)));
DeltaIL = Ij - Lj;
SR(np,loop) = norm(DeltaIL,2)/norm(Lj,2); % ||Ij-Lj||2 / ||Lj||2
end
end
save(filename,'SR','-append')
%% RMSE
for loop = 1:3
for np = 1:20
if loop == 1
    b = squeeze(abs(lcSSFP(:,:,np)));
%     b0 = squeeze(abs(lcSSFP0(:,:,np)));
elseif loop == 2
    b = squeeze(abs(idSSFP(:,:,np)));
%     b0 = squeeze(abs(idSSFP0(:,:,np)));
else
    b = squeeze(abs(ipSSFP(:,:,np)));
%     b0 = squeeze(abs(ipSSFP0(:,:,np)));
end
b = b./mean(mean(b.*mask3));%.*(imresize(phantom,192/256)~=0);
% b0 = b0./mean(mean(b0.*mask3));
% max(b(:))

% ref
a = squeeze(abs(Refimg(:,:,np)));
a = a./mean(mean(a.*mask3));%.*(imresize(phantom,192/256)~=0);

% ksp no noise
% a = b0;

RMSE(np,loop) = sqrt(mean((a(:)-b(:)).^2));
end
% figure,imshow(brighten(rot90(abs(a-b),2),0.5),[0,0.1])
figure,imshow(brighten(rot90(abs(a-b),2),0.5).*(imresize(phantom,192/256)~=0),[0,0.1]) 
% figure,imshow(brighten(rot90(abs(a-b),2),0.5),[0,0.08])
% figure,imshow(brighten(rot90(abs(a-b),2),0).*(imresize(phantom,192/256)==0),[0,0.005])
end
mean(RMSE,1)
save(filename,'RMSE','-append')
%% SNR
for loop = 1:3
for np = 1:20
if loop == 1
    a = squeeze(abs(lcSSFP(:,:,np)));
%     a0 = squeeze(abs(lcSSFP0(:,:,np)));
elseif loop == 2
    a = squeeze(abs(idSSFP(:,:,np)));
%     a0 = squeeze(abs(idSSFP0(:,:,np)));
else
    a = squeeze(abs(ipSSFP(:,:,np)));
%     a0 = squeeze(abs(ipSSFP0(:,:,np)));
end
a = a./mean(mean(a));
% a0 = a0./mean(mean(a0));

% using ref img
b = squeeze(abs(Refimg(:,:,np)));
b  = b./mean(mean(b));

% using ksp with no noise
% b = a0;

% figure,imshow(abs(a), []);
% roi = images.roi.AssistedFreehand;
% draw(roi);
% mask1 = createMask(roi);
% roi = images.roi.AssistedFreehand;
% draw(roi);
% mask2 = createMask(roi);
% mask2 = -mask2+1;
%save('F:\CVCoding\gridding_lcssfp\Simulation\ksptest48.mat', 'mask2','-append')

imgs = abs(a.*mask1);
signal  = sum(imgs(:))/sum(mask1(:)~=0);
mask2 = (imresize(phantom,192/256)==0);
a = a  - b;
imgn = abs(a.*(mask2));
% imgn = abs(a);
imgn = imgn((imgn~=0));
noise = std(imgn(:));
% SNR(np,loop)=10*log10(signal/noise);
SNR(np,loop)=(signal/noise);
end
end
save(filename,'SNR','-append')
end
%% plot UNDERSAMPLING
clear SRc SNRc RMSEc
% cd('F:\CVCoding\gridding_lcssfp\Simulation\Projection_GN90')
% cd('F:\CVCoding\gridding_lcssfp\Simulation\Gaussian Noise90')
load('F:\CVCoding\gridding_lcssfp\Simulation\Gaussian Noise90\ksptest60.mat', 'RefSNR','RefSR')
RefSNR = 1;
load('ksptest24.mat', 'SR','RMSE','SNR')
SRc(1,:,:) = SR;RMSEc(1,:,:) = RMSE;SNRc(1,:,:) = SNR./RefSNR;
load('ksptest48.mat', 'SR','RMSE','SNR')
SRc(2,:,:) = SR;RMSEc(2,:,:) = RMSE;SNRc(2,:,:) = SNR./RefSNR;
load('ksptest60.mat', 'SR','RMSE','SNR')
SRc(3,:,:) = SR;RMSEc(3,:,:) = RMSE;SNRc(3,:,:) = SNR./RefSNR;
load('ksptest80.mat', 'SR','RMSE','SNR')
SRc(4,:,:) = SR;RMSEc(4,:,:) = RMSE;SNRc(4,:,:) = SNR./RefSNR;
load('ksptest96.mat', 'SR','RMSE','SNR')
SRc(5,:,:) = SR;RMSEc(5,:,:) = RMSE;SNRc(5,:,:) = SNR./RefSNR;
load('ksptest128.mat', 'SR','RMSE','SNR')
SRc(6,:,:) = SR;RMSEc(6,:,:) = RMSE;SNRc(6,:,:) = SNR./RefSNR;
load('ksptest160.mat', 'SR','RMSE','SNR')
SRc(7,:,:) = SR;RMSEc(7,:,:) = RMSE;SNRc(7,:,:) = SNR./RefSNR;
load('ksptest192.mat', 'SR','RMSE','SNR')
SRc(8,:,:) = SR;RMSEc(8,:,:) = RMSE;SNRc(8,:,:) = SNR./RefSNR;
npj = [24,48,60,80,96,128,160,192];
npj2 = [48,60,80,96,128,160,192];
npro = transpose(sqrt(npj));
npro = ones(size(npro));
% figure,plot(npj,SRc(:,1)),hold on; plot(npj,SRc(:,2));plot(npj,SRc(:,3));title('SR')
% figure,plot(npj,RMSEc(:,1)),hold on; plot(npj,RMSEc(:,2));plot(npj,RMSEc(:,3));title('RMSE')
% figure,plot(npj,SNRc(:,1)),hold on; plot(npj,SNRc(:,2));plot(npj,SNRc(:,3));title('SNR')

figure,plot(npj,mean(SRc(:,:,1),2),'-*','LineWidth',2),
hold on; plot(npj,mean(SRc(:,:,2),2),'-*','LineWidth',2);plot(npj,mean(SRc(:,:,3),2),'-*','LineWidth',2);title('SR')
legend('lcSSFP','idSSFP','ipSSFP')
figure,plot(npj,mean(RMSEc(:,:,1),2),'-*','LineWidth',2),
hold on; plot(npj,mean(RMSEc(:,:,2),2),'-*','LineWidth',2);plot(npj,mean(RMSEc(:,:,3),2),'-*','LineWidth',2);title('RMSE')
legend('lcSSFP','idSSFP','ipSSFP'),ylim([0.25,0.8]),grid on
figure,plot(npj,mean(SNRc(:,:,1),2)./npro,'-*','LineWidth',2),
hold on; plot(npj,mean(SNRc(:,:,2),2)./npro,'-*','LineWidth',2);plot(npj,mean(SNRc(:,:,3),2)./npro,'-*','LineWidth',2);title('SNR')
legend('lcSSFP','idSSFP','ipSSFP')


figure,plot(npj2,mean(RMSEc(2:end,:,1),2),'-*','LineWidth',2),
hold on; plot(npj2,mean(RMSEc(2:end,:,2),2),'-*','LineWidth',2);plot(npj2,mean(RMSEc(2:end,:,3),2),'-*','LineWidth',2);title('RMSE')
legend('lcSSFP','idSSFP','ipSSFP')



figure,plot(npj,10.*log10(mean(SNRc(:,:,1),2)),'-*','LineWidth',2),
hold on; plot(npj,10.*log10(mean(SNRc(:,:,2),2)),'-*','LineWidth',2);plot(npj,10.*log10(mean(SNRc(:,:,3),2)),'-*','LineWidth',2);title('SNR');legend('lcSSFP','idSSFP','ipSSFP')

figure,plot(npj,mean(SNRc(:,:,1),2),'-*','LineWidth',2),
hold on; plot(npj,mean(SNRc(:,:,2),2),'-*','LineWidth',2);plot(npj,mean(SNRc(:,:,3),2),'-*','LineWidth',2);title('SNR');legend('lcSSFP','idSSFP','ipSSFP')

% % figure,plot(npj,mean(SRc(:,:,1),2)-mean(SRc(:,:,2),2)),hold on; plot(npj,mean(SRc(:,:,2),2)-mean(SRc(:,:,3),2));title('SR')
% figure,plot(npj,(mean(RMSEc(:,:,1),2)-mean(RMSEc(:,:,2),2))./mean(RMSEc(:,:,1),2),'-*','LineWidth',2),hold on; plot(npj,(mean(RMSEc(:,:,2),2)-mean(RMSEc(:,:,3),2))./(mean(RMSEc(:,:,2),2)),'-*','LineWidth',2);title('RMSE')
% legend('lcSSFP-idSSFP','idSSFP-ipSSFP')
% % figure,plot(npj,mean(SNRc(:,:,1),2)-mean(SNRc(:,:,2),2)),hold on; plot(npj,mean(SNRc(:,:,2),2)-mean(SNRc(:,:,3),2));title('SNR')
% 
% figure,plot(npj,(mean(RMSEc(:,:,1),2)-mean(RMSEc(:,:,2),2))./(mean(RMSEc(:,:,1),2)),'-*','LineWidth',2),hold on; plot(npj,(mean(RMSEc(:,:,1),2)-mean(RMSEc(:,:,3),2))./(mean(RMSEc(:,:,1),2)),'-*','LineWidth',2);
% title('RMSE');legend('lcSSFP-idSSFP','lcSSFP-ipSSFP')
