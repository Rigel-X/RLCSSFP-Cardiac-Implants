close all
clear

FinalCom = 0;
ReconCom = 1;

disp('------------------- loading the data to compare -------------------------');
% Cartesian trufi 4ch images, 30 phases, end_diastole 5 - end_systole 14 - 21
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\Metal_Carte.mat')
Metal_Carte = fliplr(Metal_Carte);

% one set of data, radial lc-ssfp, 16 phases, 3 - 7 - 13, 65 proj
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx053_10phase_90_4dyn_inter.mat')
MetalOrigin53 = four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx053_16phase_90_4dyn_inter_crop.mat')
MetalCrop53 = four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx053_16phase_90_4dyn_inter_crop_coilcut.mat')
MetalCrop_coilcut53 = four2;

% another set of data, radial lc-ssfp, 25 phases, 4 - 11 - 18, 45 proj
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx047_16phase_90_4dyn_inter_crop.mat')
MetalCrop47 = four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx047_16phase_90_4dyn_inter_crop_coilcut.mat')
MetalCrop_coilcut47 = four2;
clear four2;

if FinalCom == 1
    % compare the final results, get rid of some metal artifacts
    subplot(2,3,1); imagesc(squeeze(Metal_Carte(:,:,8)));title('Cartesian trufi','fontsize',22);colormap('gray');axis image; axis off
    subplot(2,3,2); imagesc(squeeze(MetalCrop_coilcut53(:,:,4)));title('radial lc-ssfp','fontsize',22);colormap('gray');axis image; axis off
    subplot(2,3,3); imagesc(squeeze(MetalCrop_coilcut47(:,:,6)));title('radial lc-ssfp','fontsize',22);colormap('gray');axis image; axis off
    % zoom in to show the details
    subplot(2,3,4); imagesc(squeeze(Metal_Carte(65:192,65:192,8)));colormap('gray');axis image; axis off
    subplot(2,3,5); imagesc(squeeze(MetalCrop_coilcut53(97:288,97:288,4)));colormap('gray');axis image; axis off
    subplot(2,3,6); imagesc(squeeze(MetalCrop_coilcut47(97:288,97:288,6)));colormap('gray');axis image; axis off
end

if ReconCom == 1
    % compare the different recon processing
    figure,
    subplot(1,3,1); imagesc(squeeze(MetalOrigin53(97:288,97:288,5)));title('Original Image','fontsize',22);colormap('gray');axis image; axis off
    subplot(1,3,2); imagesc(squeeze(MetalCrop53(:,:,5)));title('Cropped and Zero-filled','fontsize',22);colormap('gray');axis image; axis off
    subplot(1,3,3); imagesc(squeeze(MetalCrop_coilcut53(:,:,5)));title('Selected Coils','fontsize',22);colormap('gray');axis image; axis off
end



load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx214_10phase_40views.mat')
NotInter=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx212_10phase_40views.mat')
InterDyn=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx213_10phase_40views.mat')
InterPha=four2;
figure,
subplot(1,3,1); imshow(squeeze(abs(NotInter(30:350,30:350,8))),[0,0.1]);title('not interleaved','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,2); imshow(squeeze(abs(InterDyn(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,3); imshow(squeeze(abs(InterPha(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off



load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx220_13phase_48views.mat')
NotInter=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx216_13phase_48views.mat')
InterDyn=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx217_13phase_48views.mat')
InterPha=four2;
figure,
subplot(1,3,1); imshow(squeeze(abs(NotInter(30:350,30:350,8))),[0,0.1]);title('not interleaved','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,2); imshow(squeeze(abs(InterDyn(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,3); imshow(squeeze(abs(InterPha(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off

load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx219_30phase_49views.mat')
NotInter=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx215_30phase_49views.mat')
InterDyn=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\jx218_30phase_49views.mat')
InterPha=four2;
figure,
subplot(1,3,1); imshow(squeeze(abs(NotInter(30:350,30:350,8))),[0,0.1]);title('not interleaved','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,2); imshow(squeeze(abs(InterDyn(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics','fontsize',22);colormap('gray');axis image; axis off
subplot(1,3,3); imshow(squeeze(abs(InterPha(30:350,30:350,8))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off


load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\SamCarSAX012.mat') % 9 - 20 -
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\Sam047_17phase_60views.mat')
im47=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\Sam048_15phase_60views.mat')
im48=four2;
load('C:\Users\Jie\Desktop\CVCoding\gridding_lcssfp\Recon_Image\Sam049_22phase_60views.mat')
im49=four2;  % 5 - 10 - 
im50=squeeze(sqrt(sum(abs(imone).^2,2)));
figure,
subplot(2,3,1); imshow(squeeze(abs(a(:,:,20))),[]);title('Cartesian','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,2); imshow(squeeze(abs(im50(10,:,:))),[0,0.1]);title('one dynamic','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,4); imshow(squeeze(abs(im49(:,:,10))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,5); imshow(squeeze(abs(im48(:,:,10))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,6); imshow(squeeze(abs(im47(:,:,10))),[0,0.1]);title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off

figure,
subplot(2,3,1); imagesc(squeeze(abs(a(:,:,20))));title('Cartesian','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,2); imagesc(squeeze(abs(im50(10,:,:))));title('one dynamic','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,4); imagesc(squeeze(abs(im49(:,:,10))));title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,5); imagesc(squeeze(abs(im48(:,:,10))));title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off
subplot(2,3,6); imagesc(squeeze(abs(im47(:,:,10))));title('interleaved dynamics and phases','fontsize',22);colormap('gray');axis image; axis off

%% Compare acc and sos
% net = denoisingNetwork('DnCNN');
% denoisedI = denoiseImage(noisyI,net);
% thePSNR = psnr(testImage, referenceImage);
% theMSE = immse(testImage, referenceImage);
for nphase = 1:phases
    noisyI = squeeze(abs(four1(:,:,nphase)));
    Ref1(:,:,nphase) = denoiseImage(noisyI,net);
    PSNR1(nphase) = psnr(noisyI, squeeze(Ref1(:,:,nphase)));
    MSE1(nphase) = immse(noisyI, squeeze(Ref1(:,:,nphase)));
    noisyI = squeeze(abs(four2(:,:,nphase)));
    Ref2(:,:,nphase) = denoiseImage(noisyI,net);
    PSNR2(nphase) = psnr(noisyI, squeeze(Ref2(:,:,nphase)));
    MSE2(nphase) = immse(noisyI, squeeze(Ref2(:,:,nphase)));    
end
for nphase = 1:phases
    Ij = squeeze(four2(:,:,nphase));
    IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
    Lj = fftshift(fft2(ifftshift(IjTemp)));
    DeltaIL = Ij - Lj;
    SR2(nphase) = norm(DeltaIL,2)/norm(Lj,2); % ||Ij-Lj||2 / ||Lj||2
end
for nphase = 1:phases
    Ij = squeeze(four1(:,:,nphase));
    IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
    Lj = fftshift(fft2(ifftshift(IjTemp)));
    DeltaIL = Ij - Lj;
    SR1(nphase) = norm(DeltaIL,2)/norm(Lj,2); % ||Ij-Lj||2 / ||Lj||2
end
len=size(four2,3);
for j = 1:len
    % CV gif
    subplot(2,3,1)
    sr=SR1(j);psnr=PSNR1(j);mse=MSE1(j);
    imshow(rot90(abs(four1(:,:,j))),[0,0.15]);title(['sos phase=',int2str(j),', SR=',num2str(sr),', PSNR=',num2str(psnr),', MSE=',num2str(mse)]);
    axis square,colormap('gray'),colorbar;
    subplot(2,3,2)
    sr=SR2(j);psnr=PSNR2(j);mse=MSE2(j);
    imshow(rot90(abs(four2(:,:,j))),[0,0.15]);title(['acc phase=',int2str(j),', SR=',num2str(sr),', PSNR=',num2str(psnr),', MSE=',num2str(mse)]);
    axis square,colormap('gray'),colorbar;
    subplot(2,3,4)
    imshow(rot90(abs(Ref1(:,:,j)))-rot90(abs(four1(:,:,j))),[]);title(['denoised reference image, sos phase=',int2str(j)]);
    axis square,colormap('gray'),colorbar;
    subplot(2,3,5)
    imshow(rot90(abs(Ref2(:,:,j)))-rot90(abs(four2(:,:,j))),[]);title(['denoised reference image, acc phase=',int2str(j)]);
    axis square,colormap('gray'),colorbar;
    subplot(2,3,6)
    imshow(rot90(abs(four1(:,:,j))-abs(four2(:,:,j))),[]);title(['sos-acc phase=',int2str(j)]);
    axis square,colormap('gray'),colorbar;
    set(gcf,'Color',[1 1 1]);
    set(gcf, 'Position', get(0, 'Screensize'));
    nn=getframe(gcf);
    im=frame2im(nn);
    [I,map]=rgb2ind(im,256);
    if j==1
        imwrite(I,map,'Recon_Image\anj\Anj058__.gif','gif','loopcount',inf,'Delaytime',2)
    else
        imwrite(I,map,'Recon_Image\anj\Anj058__.gif','gif','writemode','append','Delaytime',2)
    end
end



%% Compare with DC denoising
for nphase = 1:phases
    noisyI = squeeze(abs(four2(:,:,nphase)));
    Ref2(:,:,nphase) = denoiseImage(noisyI,net);
    MSE(nphase) = immse(noisyI, squeeze(Ref2(:,:,nphase)));
end
len=size(four2,3);
for j = 1:len
    % CV gif
    subplot(1,3,1)
    imshow(abs(four2(:,:,j)),[0,0.05]);title(['original phase=',int2str(j)]);
    axis square,colormap('gray'),colorbar;
    subplot(1,3,2)
    imshow(abs(Ref2(:,:,j)),[0,0.05]);title(['DC denoised phase=',int2str(j)]);
    axis square,colormap('gray'),colorbar;
    subplot(1,3,3)
    mse=MSE(j);
    imshow(abs(Ref2(:,:,j))-abs(four2(:,:,j)),[]);title(['denoise - original phase=',int2str(j),', MSE=',num2str(mse)]);
    axis square,colormap('gray'),colorbar;
    set(gcf,'Color',[1 1 1]);
    set(gcf, 'Position', get(0, 'Screensize'));
    nn=getframe(gcf);
    im=frame2im(nn);
    [I,map]=rgb2ind(im,256);
    if j==1
        imwrite(I,map,'Recon_Image\rongtao\rt047.gif','gif','loopcount',inf,'Delaytime',2)
    else
        imwrite(I,map,'Recon_Image\rongtao\rt047.gif','gif','writemode','append','Delaytime',2)
    end
end


%% Compare different recon methods
%load('Sam049nufft.mat')
four2=img_nufft1+img_nufft2+img_nufft3+img_nufft4;
len=size(four2,3);
for j = 1:len
imshow(abs(four2(:,:,j)),[0,0.0035]);
title(['nufft phase=',int2str(j)]);colorbar;
%imshow(abs(four2(97:288,97:288,j)),[]);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'sam047nufft.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'sam047nufft.gif','gif','writemode','append','Delaytime',0.05)
end
end

load('Sam047_nlinv.mat')
four2=img_nlinv1+img_nlinv2+img_nlinv3+img_nlinv4;
four2=rot90(four2,2);
for j = 1:len
imshow(abs(four2(:,:,j)),[0,0.6]);
title(['nlinv phase=',int2str(j)]);colorbar;
%imshow(abs(four2(97:288,97:288,j)),[]);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'sam047nlinv.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'sam047nlinv.gif','gif','writemode','append','Delaytime',0.05)
end
end

load('Sam047_pics.mat')
four2=reco_l11+reco_l12+reco_l13+reco_l14;
four2=rot90(four2,2);
for j = 1:len
imshow(abs(four2(:,:,j)),[0,0.025]);
title(['pics L1 wavelet phase=',int2str(j)]);colorbar;
%imshow(abs(four2(97:288,97:288,j)),[]);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'sam047picsL1.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'sam047picsL1.gif','gif','writemode','append','Delaytime',0.05)
end
end

four2=reco_l21+reco_l22+reco_l23+reco_l24;
four2=rot90(four2,2);
for j = 1:len
imshow(abs(four2(:,:,j)),[0,50]);
title(['pics L2 phase=',int2str(j)]);colorbar;
%imshow(abs(four2(97:288,97:288,j)),[]);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'sam047picsL2.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'sam047picsL2.gif','gif','writemode','append','Delaytime',0.05)
end
end

four2=reco_tv1+reco_tv2+reco_tv3+reco_tv4;
four2=rot90(four2,2);
for j = 1:len
% CV gif
imshow(abs(four2(:,:,j)),[]);
title(['pics total variation phase=',int2str(j)]);colorbar;
%imshow(abs(four2(97:288,97:288,j)),[]);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'sam047picstv.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'sam047picstv.gif','gif','writemode','append','Delaytime',0.05)
end
end






