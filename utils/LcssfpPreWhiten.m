% 2022/3/17, PreWhitening of kspace
% jie xiang @yale mrrc

function [ksp_prewhitened,dataMnew] = LcssfpPreWhiten(NoiseData,coils,ksp_raw)
%JX_PREWHITEN Summary of this function goes here
%   prewhitening before recon, in case of broken coils
%   ksp_raw: original kspace
%   coils: number of coils
%   NoiseData: turning down the RF for noise cov matrix calculation

NoiseData = NoiseData(1:coils*128,:);
for nc = 1:coils
for nl = 1:128 % lines
noisedata(nl,:,nc) = NoiseData((nl-1)*coils+nc,:);
end
end
eta = permute(noisedata,[3,1,2]);
eta = reshape(eta,[nc,nl*512]);
Nsamples = size(eta,2);
psi = (1/(Nsamples-1))*(eta*eta');
% figure,imshow(abs(psi),[]),colormap('jet')

%%
L = chol(psi,'lower');
L_inv = inv(L);
etanew = L_inv * eta;
dataMnew  = (1/(Nsamples-1))*(etanew*etanew');
figure,
subplot(1,2,1),imshow(abs(psi),[]),colormap('jet'),title('original psi')
subplot(1,2,2),imshow(abs(dataMnew),[]),colormap('jet'),title('pre-whitened psi')

[nc,nphase,nPE,nFE] = size(ksp_raw);
for np = 1:nphase
    data0 = squeeze(ksp_raw(:,np,:,:));
% figure,imshow(abs(squeeze(data(1,:,:))),[])
data = reshape(data0,[nc,nPE*nFE]);
% dataM1  = (1/(Nsamples-1))*(data*data');
% figure,imshow(abs(dataM1),[]),colormap('jet')
%%
data = L_inv * data;
% dataM2  = (1/(Nsamples-1))*(data*data');
% figure,
% subplot(1,2,1),imshow(abs(dataM1),[]),colormap('jet'),title('original data')
% subplot(1,2,2),imshow(abs(dataM2),[]),colormap('jet'),title('pre-whitened data')
% csm = L_inv * csm; % coil sensitivity map
%%
ksp_prewhitened = zeros(nc,np,nPE,nFE);
ksp_prewhitened(:,np,:,:) = reshape(data,[nc,nPE,nFE]);
end
% max0 = max(abs(data0(:)));
% max1 = max(abs(data(:)));
% figure,
% subplot(3,1,1),imshow(squeeze(abs(ksp_raw(1,np,:,:)))./max0,[]),title('original ksp')
% subplot(3,1,2),imshow(squeeze(abs(ksp_prewhitened(1,np,:,:)))./max1,[]),title('pre-whitened ksp')
% subplot(3,1,3),imshow(100*(squeeze(abs(ksp_prewhitened(1,np,:,:)))./max1-squeeze(abs(ksp_raw(1,np,:,:)))./max0)./(squeeze(abs(ksp_raw(1,np,:,:)))./max0),[]),title('ksp difference %')
end

