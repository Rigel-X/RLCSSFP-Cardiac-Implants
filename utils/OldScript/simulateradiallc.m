%SIMULATE ARTEFACTS FROM RADIAL LCSSFP

%SELECT A ROI FROM AN IMAGE OR SIMULATED IMAGE
% figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,1,:,:))));
% colormap(gray); imagesc(flipud((finalIm(:,:)))); axis equal
% ph1=roipoly();
% ph1=double(ph1);
% [Nx,Ny]=size(ph1);
% measures = 4; %number of passes

%OR LOAD ROI
% load ph1; 
%LOAD SIGNAL VS FREQUENCY DATA OR GENERATE
load store;

[Nx,Ny]=size(ph1);
freqoff=repmat(linspace(-500,500,Ny),Nx,1); %freq map that goes from 500 to -500. 
% figure;imagesc(freqoff)

%CALCULATE SIGNAL IN PIXEL CORRESPONDING TO FREQUENCY
for p = 1:measures
for i=1:Nx,
    for j=1:Ny,
        freq=ceil(freqoff(i,j))+501;
        signal = sqrt((store(p,1,freq)+(store(p,2,freq))).^2);
        si = store(p,1,freq)+sqrt(-1)*store(p,2,freq);
        phase = angle(si);
        phantom(p,i,j)=ph1(i,j) *signal*(cos(phase) +sqrt(-1)*sin(phase)); %*mag(freq)*phase(freq);  % the signal and phase of ph1 w
    end
end

Np=250; %NUMBER OF UNDERSAMPLED PROJECTIONS
off=0;
%SIMULATE INTERLEAVING COMMENT OUT FOR NO-INTERLEAVING
if p == 1; q = 0 
else if p == 2; q = 1
    else if p == 3; q =2
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
imager=iradon(radonr,angles,Nx,Ny);
imagei=iradon(radoni,angles,Nx,Ny);
newim=complex(imager,imagei);

allims(p,:,:) = newim;

% figure, imagesc(squeeze(abs(radonx)));
% title('fbp of data phase cycling 1');
% figure, imagesc(squeeze((imager)));
% figure, imagesc(squeeze((imagei)));
% figure, imagesc(angle((newim)));
% title(sprintf( "%s%d", 'Phase Cyling = ', p))
end

%COMBINE IMAGES FROM EACH PASS OF SSFP
combinationim = squeeze(allims(1,:,:) + allims(2,:,:) + allims(3,:,:) + allims(4,:,:));
scalled = combinationim/max(max(combinationim)); %SCALE FROM 0 TO 1

% figure, colormap(gray); imagesc(abs(scalled))
% title(sprintf("%s%s%d%s%d", 'Combination - Interleaved std Order', "p = ", p, "Np = ", Np))
% figure, colormap(gray);imagesc(abs(squeeze(allims(3,:,:))))
% title('0-180')
% figure, colormap(gray); imagesc(abs(squeeze(allims(1,:,:))))
% title('0-0')


%SELECT A SIGNAL ROI OR LOAD 
%CALCULATE PROPERTIES
%     siggy = roipoly;
%     signal_filter = double(siggy);
% load signal_filter
    signal_filter(signal_filter == 0) = NaN;
    shnigal = scalled.*signal_filter;
    sig_mean1 = mean(shnigal(~isnan(shnigal)));
    sig_max = max(shnigal(~isnan(shnigal)));
    sig_min = min(shnigal(~isnan(shnigal)));
    sig_total = sum(shnigal(~isnan(shnigal)));
    sig_sd = std(shnigal(~isnan(shnigal)));
%     save signal_filter signal_filter

%SELECT OR LOAD NOISE/ARTEFACT ROI
%CALCULATE PROPERTIES
%     noise = roipoly;
% load noise_filter
%     noise_filter = double(noise);
    noise_filter(noise_filter == 0) = NaN;
    naise = scalled.*noise_filter;
    noi_mean1 = mean(naise(~isnan(naise)));
    noi_max = max(naise(~isnan(naise)));
    noi_min = min(naise(~isnan(naise)));
    noi_total = sum(naise(~isnan(naise)));
    noi_sd = std(naise(~isnan(naise)));
    
%     save noise_filter noise_filter

%DISPLAY PROPERTIES
words = "The max signal is: %0.5f  \nThe min signal is: %0.5f \nThe mean signal  is: %0.5f \nThe standard deviation of signal is: %0.5f \nThe total signal is: %0.5f ";
words2 ="The max artefact is: %0.5f \nThe min artefact is: %0.5f \nThe mean artefact  is: %0.5f  \nThe standard deviation of artefact is: %0.5f \nThe total artefact is: %0.5f";
sprintf(words,abs(sig_max),abs(sig_min),abs(sig_mean1),abs(sig_sd),abs(sig_total))
sprintf(words2, abs(noi_max), abs(noi_min),abs(noi_mean1), abs(noi_sd), abs(noi_total))