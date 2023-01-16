%SIMULATE ARTEFACTS FROM RADIAL LCSSFP

%SELECT A ROI FROM AN IMAGE OR SIMULATED IMAGE
% figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,1,:,:))));
% colormap(gray); imagesc((fliplr(finalIm(:,:)))); axis equal
% ph1=roipoly();
% ph1=double(ph1);
% [Nx,Ny]=size(ph1);
measures = 4; %number of passes


for disSiz = 192%[160,192,256]
close all;
clear arr2
load simulationRois.mat
Npidx=1;
%OR LOAD ROI
load ph1; 
%LOAD SIGNAL VS FREQUENCY DATA OR GENERATE
% dim1: pass; dim2: real and imag parts, dim3: 1001 frequencies (from
% -500hz to 500hz)
load storeT1_1600_T2_200_fa_30.mat;

load noiseForSimulationNew.mat
ph1 = ph1(33:147,56:204);
ph1 = padarray(ph1,[10,0]);
% ph1 = interp2(ph1,1);
[xx,yy]=meshgrid(linspace(0,1,size(ph1,2)),linspace(0,1,size(ph1,1)));
% xx = mat2gray(xx); yy = mat2gray(yy);
[xx2,yy2] = meshgrid(linspace(0,1,192));
% xx2 = mat2gray(xx2); yy2 = mat2gray(yy2);
ph1 = interp2(xx,yy,ph1,xx2,yy2);
% ph1 = ph1(round(size(ph1,1)/2)-96:round(size(ph1,1)/2)+95,...
%     round(size(ph1,2)/2)-96:round(size(ph1,2)/2)+95);
xsmall = -96:95;
[X1, Y1] = meshgrid(xsmall/192);
[X2, Y2] = meshgrid((-disSiz/2:(disSiz/2-1))/disSiz);
ph1 = interp2(X1, Y1, ph1, X2, Y2, 'linear', 0);
% ph1 = reshape(ph1, [disSiz, disSiz]);
figure; imagesc(ph1);
[Nx,Ny]=size(ph1);
numpix = length(find(ph1>0.02));

freqoff=repmat(linspace(-1500,1500,Ny),Nx,1); %freq map that goes from 500 to -500. 
figure;imagesc(freqoff);
figure; imagesc(ph1); axis equal off
figure; imagesc(freqoff.*double(ph1>0)); axis equal off

lcArr = cell(1,1,12);
lcMsArr = lcArr;
lcSosArr = lcArr;
%
for Np=[16,32,48,64,96,128,160,192,256,320,384]
    %%
    fprintf('%d: iter\n', Np);
for noiselevel = 2

%%
%CALCULATE SIGNAL IN PIXEL CORRESPONDING TO FREQUENCY
phantom = zeros(measures, Nx, Ny);
phantomNoise = phantom;
allims = zeros(measures, Ny, Ny);
allims2 = allims;
for p = 1:measures
    for i=1:Nx,
        for j=1:Ny,
            freq=ceil(mod(freqoff(i,j)+500,1000)-500)+501;
%             signal = sqrt((store(p,1,freq)+store(p,2,freq)).^2);
%             si = store(p,1,freq)+1i*store(p,2,freq);
%             phase = angle(si);
%             phantom(p,i,j)=ph1(i,j) *signal*(cos(phase) +1i*sin(phase)); %*mag(freq)*phase(freq);  % the signal and phase of ph1 w
            phantom(p,i,j) = ph1(i,j) * (store(p,1,freq) + 1i * store(p,2,freq));
        end
    end
    
    figure; subplot(121);imagesc(abs(squeeze(phantom(p,:,:))));
    
%     phantom(p,:,:) = phantom(p,:,:) + noiseForSimulationNew(p,1:Nx,1:Nx)/noiselevel;

     %NUMBER OF UNDERSAMPLED PROJECTIONS
%     off=0;
    %SIMULATE INTERLEAVING COMMENT OUT FOR NO-INTERLEAVING
    if p == 1; q = 0 
    else if p == 2; q = 2
        else if p == 3; q =1
            else if p == 4; q = 3
    
                end
            end
        end
    end
    %%%%%
    %GENERATE PROJECTIONS ANGLES AND RECONSTRUCT DATA RADIALLY
    off=(180/(Np*4))*q; %offset angle can be (180/(Np*4))*1 or 2 or 3
    angles=[off:180/Np:180-180/Np+off];
    
    radonr=radon(real(squeeze(phantom(p,:,:))),angles);
    radoni=radon(imag(squeeze(phantom(p,:,:))),angles);
    
    rnorm = mean(abs(radonr(:)+1i*radoni(:)));
    nstd = rnorm/1/noiselevel/sqrt(2)*0;
    noise = nstd*randn(size(radonr)) + 1i * nstd *randn(size(radonr));
    radonr = radonr + real(noise);
    radoni = radoni + imag(noise);
    
    
    imager=iradon(radonr,angles,Nx,Ny);
    imagei=iradon(radoni,angles,Nx,Ny);
    newim=complex(imager,imagei);

    allims(p,:,:) = newim;
    subplot(122); imagesc(abs(squeeze(allims(p,:,:))));
    
    
    % no interleaving
    off=0; %offset angle can be (180/(Np*4))*1 or 2 or 3
    angles=[off:180/Np:180-180/Np+off];
    
    radonr=radon(real(squeeze(phantom(p,:,:))),angles);
    radoni=radon(imag(squeeze(phantom(p,:,:))),angles);
    
    radonr = radonr + real(noise);
    radoni = radoni + imag(noise);
    
    
    imager=iradon(radonr,angles,Nx,Ny);
    imagei=iradon(radoni,angles,Nx,Ny);
    newim=complex(imager,imagei);

    allims2(p,:,:) = newim;
    % figure, imagesc(squeeze(abs(radonx)));
    % title('fbp of data phase cycling 1');
    % figure, imagesc(squeeze((imager)));
    % figure, imagesc(squeeze((imagei)));
    % figure, imagesc(angle((newim)));
    % title(sprintf( "%s%d", 'Phase Cyling = ', p))
end
%%
%COMBINE IMAGES FROM EACH PASS OF SSFP
lc = squeeze(allims(1,:,:) + allims(2,:,:) + allims(3,:,:) + allims(4,:,:));
lcMs = squeeze(abs(allims(1,:,:)) + abs(allims(2,:,:)) + abs(allims(3,:,:)) + abs(allims(4,:,:)));
lcSos = sqrt(squeeze(abs(allims(1,:,:)).^2 + abs(allims(2,:,:)).^2 + abs(allims(3,:,:)).^2 + abs(allims(4,:,:)).^2));
lc2 = squeeze(allims2(1,:,:) + allims2(2,:,:) + allims2(3,:,:) + allims2(4,:,:));
lc2Ms = squeeze(abs(allims2(1,:,:)) + abs(allims2(2,:,:)) + abs(allims2(3,:,:)) + abs(allims2(4,:,:)));
lc2Sos = sqrt(squeeze(abs(allims2(1,:,:)).^2 + abs(allims2(2,:,:)).^2 + abs(allims2(3,:,:)).^2 + abs(allims2(4,:,:)).^2));
% scalled = combinationim/max(max(combinationim)); %SCALE FROM 0 TO 1

figure, colormap(gray); imagesc(abs(lc))

%% measure SNR, banding artifact, and streak artifact
if 0%Np==16
    figure; imagesc(abs(squeeze(allims(3,:,:))));
    obj = roipoly;
    objcenter = roipoly;
    bk = roipoly;
    save('simulationROIs.mat', 'obj', 'objcenter', 'bk');
end
load simulationROIs.mat
obj = interp2(X1, Y1, double(obj), X2, Y2, 'linear', 0);
obj = logical(round(obj)==1);
objcenter = interp2(X1, Y1, double(objcenter), X2, Y2, 'linear', 0);
objcenter = logical(round(objcenter));
bk = interp2(X1, Y1, double(bk), X2, Y2, 'linear', 0);
bk = logical(round(bk));
% angles=[0:180/192:180-180/192];
% radonr=radon(real(ssfp),angles);
% radoni=radon(imag(ssfp),angles);
% imager=iradon(radonr,angles,Nx,Ny);
% imagei=iradon(radoni,angles,Nx,Ny);
% ssfp=abs(complex(imager,imagei));
lc = abs(lc);
lcMs = abs(lcMs);
lcSos = abs(lcSos);
lc2 = abs(lc2);
SigMean = [mean(lc(obj)), mean(lcSos(obj)), mean(lcMs(obj)), mean(lc2(obj)), mean(lc2Sos(obj)), mean(lc2Ms(obj))];
SigMean2 = [mean(lc(objcenter)), mean(lcSos(objcenter)), mean(lcMs(objcenter)), mean(lc2(objcenter)), mean(lc2Sos(objcenter)), mean(lc2Ms(objcenter))];
SigStd = [std(lc(obj)), std(lcSos(obj)), std(lcMs(obj)), std(lc2(obj)), std(lc2Sos(obj)), std(lc2Ms(obj))];
NoiseStd = [std(lc(bk)), std(lcSos(bk)), std(lcMs(bk)), std(lc2(bk)), std(lc2Sos(bk)), std(lc2Ms(bk))];
NoiseStd2 = [norm(lc(bk)), norm(lcSos(bk)), norm(lcMs(bk)), norm(lc2(bk)), norm(lc2Sos(bk)), norm(lc2Ms(bk))];
NoiseStd2 = NoiseStd2/numpix;
SNR = SigMean2./NoiseStd;
SNR2 = SigMean2./NoiseStd2;
SNR2 = SNR2(:);

lcArr{Npidx} = lc;
lcMsArr{Npidx} = lcMs;
lcSosArr{Npidx} = lcSos;

band = SigMean./SigStd;
SNR = SNR(:);
band = band(:);
arr = [SNR,band, SNR2];
arr2{noiselevel, Npidx} = arr;

% save(['paperSimulationLCSSFP_Np=' num2str(Np) 'NoiseLevel=' num2str(noiselevel) '.mat'], 'lc', 'lcMs', 'lcSos', 'lc2', 'lc2Ms', 'lc2Sos');

if Np==Nx
    ssfparr{noiselevel} = squeeze(allims(3,:,:));
    SNRsa(noiselevel) = mean(abs(ssfparr{noiselevel}(objcenter)))/std(abs(ssfparr{noiselevel}(bk)));
    BandSA(noiselevel) = mean(abs(ssfparr{noiselevel}(obj)))/std(abs(ssfparr{noiselevel}(obj)));
end

end
Npidx = Npidx+1;
end
arr3 = cell2mat(arr2);
lcArr = cell2mat(lcArr);
lcMsArr = cell2mat(lcMsArr);
lcSosArr = cell2mat(lcSosArr);

 %% generate graphs and images
for noiselevel = 2%1:3
close all
nparr=[16,32,48,64,96,128,160,192,256,320,384];
SNRarr = arr3(:,1:3:end);
bandarr = arr3(:,2:3:end);
SNR2arr = arr3(:,3:3:end);
% ssfp = abs(squeeze(allims(3,:,:)));
% SNRssfp = mean(ssfp(objcenter))/std(ssfp(bk));
% Bandssfp = mean(ssfp(obj))/std(ssfp(obj));

% 1. show criteria changes with Np at NoiseLevel=2
c = {'-xb','-xr','-xg','-.ob','-.or','-.og','--^k'};
cf = [0.6551,0.6551,0.6551,0.6551];
cf = [1, 1, 1, 1];
ind = 0;%(noiselevel-1)*6;
for idx = 1:6
    figure(1); plot(nparr,SNRarr(ind+idx,:),c{idx});hold on; 
%     if idx<5
%         figure(2); plot(nparr,SNRarr(0+idx,:)./sqrt(nparr*4),c{idx});hold on; 
%     else
%         figure(2); plot(nparr,SNRarr(0+idx,:)./sqrt(Nx),c{idx});hold on; 
%     end
    figure(2); plot(nparr,bandarr(ind+idx,:),c{idx});hold on; 
end
figure(1);plot(nparr, SNRsa(noiselevel),c{end});
figure(2);plot(nparr, BandSA(noiselevel)*ones(size(nparr)),c{end});
figure(1);legend('CS-Intd', 'SOS-Intd', 'MS-Intd', 'CS-Nonintd', 'SOS-Nonintd', 'MS-Nonintd', 'F.S. bSSFP'); hold off;
% figure(2);legend('CS LC-SSFP', 'MS LC-SSFP', 'SOS LC-SSFP', 'NonIntd. LC-SSFP', 'bSSFP'); hold off;
figure(2);legend('CS-Intd', 'SOS-Intd', 'MS-Intd', 'CS-Nonintd', 'SOS-Nonintd', 'MS-Nonintd', 'F.S. bSSFP'); hold off;

for idx = 1:6
    figure(3); plot(nparr,SNRarr(ind+idx,:)./SNRsa(noiselevel), c{idx}); hold on;
    figure(4); plot(nparr, (SNRarr(ind+idx,:)./sqrt(nparr*4))/(SNRsa(noiselevel)/sqrt(Nx)), c{idx}); hold on;
    figure(5); plot(nparr, bandarr(ind+idx,:)/BandSA(noiselevel), c{idx}); hold on;
end
figure(3);legend('CS-Intd / bSSFP', 'SOS-Intd / bSSFP', 'MS-Intd / bSSFP', 'CS-Nonintd / bSSFP', 'SOS-Nonintd / bSSFP', 'MS-Nonintd / bSSFP'); hold off;
figure(4);legend('CS-Intd / bSSFP', 'SOS-Intd / bSSFP', 'MS-Intd / bSSFP', 'CS-Nonintd / bSSFP', 'SOS-Nonintd / bSSFP', 'MS-Nonintd / bSSFP'); hold off;
figure(5);legend('CS-Intd / bSSFP', 'SOS-Intd / bSSFP', 'MS-Intd / bSSFP', 'CS-Nonintd / bSSFP', 'SOS-Nonintd / bSSFP', 'MS-Nonintd / bSSFP'); hold off;
%save image
% for fidx = 1:5
%     figure(fidx);savefig(['fig' num2str(fidx) '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.fig']);
%     export_fig(['fig' num2str(fidx) '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.jpg'], '-a1');
% end
% part 2 criteria over noiseStd
% c = {'-xb','-xg','-xm','-ok','-^r'};
% for idx = 1:5
%     figure(1); plot(snrlevel,SNRarr((0:5:10)+idx,3),c{idx});hold on;     
%     figure(3); plot(snrlevel,bandarr((0:5:10)+idx,3),c{idx});hold on; 
% end
% figure(1);legend('CS LC-SSFP', 'MS LC-SSFP', 'SOS LC-SSFP', 'NonInt. LC-SSFP', 'bSSFP'); hold off;
% % figure(2);legend('CS LC-SSFP', 'MS LC-SSFP', 'SOS LC-SSFP', 'NonInt. LC-SSFP', 'bSSFP'); hold off;
% figure(3);legend('CS LC-SSFP', 'MS LC-SSFP', 'SOS LC-SSFP', 'NonInt. LC-SSFP', 'bSSFP'); hold off;

% save(['simulationResult01292019_ImageSize=' num2str(Nx) '.mat']); 
end
%% part 3. show images at Np=48 for the 3 noise levels
% close all
% r = 1;
% name = {'paperSimulationLCSSFP_Np=24NoiseLevel=1',...
%     'paperSimulationLCSSFP_Np=48NoiseLevel=1'};
% 
% for idx = 1:length(name)
%     
%     load([name{idx} '.mat']);
% 
%     figure; imagesc((lc), [0,max(lc(:))*r]);colormap(gray);axis equal off; title(['cslc' name{idx}]);
%     export_fig(['cs-lc' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((lcMs), [0,max(lcMs(:))*r]);colormap(gray);axis equal off; title(['mslc' name{idx}]);
%     export_fig(['ms-lc' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((lcSos), [0,max(lcSos(:))*r]);colormap(gray);axis equal off; title(['soslc' name{idx}]);
%     export_fig(['sos-lc' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((lc2), [0,max(lc2(:))*r]);colormap(gray);axis equal off; title(['nonint' name{idx}]);
%     export_fig(['nonint' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((lc2Ms), [0,max(lc2Ms(:))*r]);colormap(gray);axis equal off; title(['nonintMs' name{idx}]);
%     export_fig(['nonint-Ms' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((lc2Sos), [0,max(lc2Sos(:))*r]);colormap(gray);axis equal off; title(['nonintSos' name{idx}]);
%     export_fig(['nonint-Sos' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
%     figure; imagesc((ssfp), [0,max(ssfp(:))*r]);colormap(gray);axis equal off; title(['bssfp' name{idx}]);
%     export_fig(['bssfp' name{idx} '_NoiseLevel' num2str(noiselevel) '_imagesize' num2str(Nx) '.tif'], '-r300', '-a1');
% 
% end


end



% save('newSimulationEverything.mat');

% title(sprintf("%s%s%d%s%d", 'Combination - Interleaved std Order', "p = ", p, "Np = ", Np))
% figure, colormap(gray);imagesc(abs(squeeze(allims(3,:,:))))
% title('0-180')
% figure, colormap(gray); imagesc(abs(squeeze(allims(1,:,:))))
% title('0-0')


% 
% %SELECT A SIGNAL ROI OR LOAD 
% %CALCULATE PROPERTIES
% %     siggy = roipoly;
% %     signal_filter = double(siggy);
% % load signal_filter
%     signal_filter(signal_filter == 0) = NaN;
%     shnigal = scalled.*signal_filter;
%     sig_mean1 = mean(shnigal(~isnan(shnigal)));
%     sig_max = max(shnigal(~isnan(shnigal)));
%     sig_min = min(shnigal(~isnan(shnigal)));
%     sig_total = sum(shnigal(~isnan(shnigal)));
%     sig_sd = std(shnigal(~isnan(shnigal)));
% %     save signal_filter signal_filter
% 
% %SELECT OR LOAD NOISE/ARTEFACT ROI
% %CALCULATE PROPERTIES
% %     noise = roipo