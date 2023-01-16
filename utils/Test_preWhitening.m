clc
close all
rawdata = MultiRAID(1).rawdata;
rawdata = rawdata(1:8192,:);
for nc = 1:coils
for nl = 1:128 % lines
noiseData(nl,:,nc) = rawdata((nl-1)*64+nc,:);
end
end
eta = permute(noiseData,[3,1,2]);
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

data0 = squeeze(first(:,1,:,:));
% figure,imshow(abs(squeeze(data(1,:,:))),[])
[nc,nPE,nFE] = size(data0);
data = reshape(data0,[nc,nPE*nFE]);
dataM1  = (1/(Nsamples-1))*(data*data');
% figure,imshow(abs(dataM1),[]),colormap('jet')
%%
data = L_inv * data;
dataM2  = (1/(Nsamples-1))*(data*data');
figure,
subplot(1,2,1),imshow(abs(dataM1),[]),colormap('jet'),title('original data')
subplot(1,2,2),imshow(abs(dataM2),[]),colormap('jet'),title('pre-whitened data')
% csm = L_inv * csm; % coil sensitivity map
%%
data = reshape(data,[nc,nPE,nFE]);
max0 = max(abs(data0(:)));
max1 = max(abs(data(:)));
figure,
subplot(3,1,1),imshow(squeeze(abs(data0(1,:,:)))./max0,[]),title('original ksp')
subplot(3,1,2),imshow(squeeze(abs(data(1,:,:)))./max1,[]),title('pre-whitened ksp')
subplot(3,1,3),imshow(100.*(squeeze(abs(data0(1,:,:)))./max0-squeeze(abs(data(1,:,:)))./max1)./(squeeze(abs(data0(1,:,:)))./max0),[]),title('ksp difference, %')