%SIMULATING DYPR
%DYNAMICALLY PHASE CYCLED RADIAL SSFP

%SPECIFY T1 T2
myoT1 = 1500; 
myoT2 = 200;  

%SPECIFY PARAMETERS
TR = 3; %4 6 10];   %ms
acq = 800;          %length of acquisition ms
T = acq;
proj = 250;         %number for (undersampled) projections i.e. Np/4
meas = 4;       
dummy = 1000;       %number of dummy scans
totp = meas*proj;   %Np
ns = totp + dummy;  %Total number of excitations 
alpha = 60;         %flip angle 
df =[-500:500];     %off resonance frequency range Hz

magnetization=zeros(length(df), 2, ns);

%INITIALIZE
mb = [];            %magnetisation blood
mm = mb;            %magnetisation myocardium
t = zeros(1,length(mb));
mb(:,1) = [0;0;1];  %set column one to 0, 0, 1, ie mz = 1
mm(:,1) = [0;0;1];  
t(1) = 0;
idx = 2;
k=1;
signal=zeros(length(df),2,ns);

%GENERATE PHASE INCREMENT APPLIED BETWEEN EACH EXCITATION
pha = 0;
for s=2:ns
pha(s)=((360/(totp))+pha(s-1));
end

for f = 1:length(df)
    numTR = ns;                 %Repeats per acquisition
    e2m = exp(-TR/myoT2);       %standard exponential terms 
    e1m = exp(-TR/myoT1);
    %CALL FUNCTION TO CALCULATE SIGNAL HEURISTICALLY 
    [magnetization] = dypr (TR, df(f), e1m, e2m, alpha, idx, numTR, mm, t, myoT1, myoT2, ns,f, pha);
    signal(f,1,:) = magnetization(f,1,:); 
    signal(f,2,:) = magnetization(f,2,:);
%     for s =1000:ns
%     xy(f,s) = (sqrt(signal(f,1,s)).^2 + (signal(f,2,s)).^2);
%     end
end
%PLOTTING
% figure;
% for s=ns/2:ns
%     hold off
% plot(df,abs(xy(:,s)));hold on
% xlabel('Frequency(Hz)')
% ylabel('Magnetization')
% title(sprintf( "%s%d%s%d", 'Spectrum for projection =  ',  s, 'Phase increment = ', ph(s)));
% pause(0.5);
% end
% for f = 1:length(df)
%     total(f) = sum(squeeze(xy(f,:)));
% end 
% figure;plot(df,abs(total))
% title('Summed signal from all projections')
% xlabel('Frequency')
% ylabel('Magnetization')
% 
% figure; plot((squeeze(signal(ceil(length(df)/2), 1,1000:ns))), squeeze((-1)*signal(ceil(length(df)/2), 2,1000:ns)))
% xlabel('Mx (Real Mxy)')
% ylabel('My (Imag Mxy)')



% SIMULATE STREAKING ARTEFACT

%CHOOSE ROI FROM GRE IMAGE
% figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,1,:,:))));
%     colormap(gray); imagesc(flipud((finalIm(:,40:186)))); axis equal
% ph1=roipoly();
% ph1=double(ph1);
% [Nx,Ny]=size(ph1);

clear combinationim  magnetization phantom angles radonr radoni imager imagei 
load ph1 
[Nx,Ny]=size(ph1);
freqoff=repmat(linspace(-500,500,Ny),Nx,1); % generate freq map that goes from 500 to -500. 

phantom = zeros(ns-1000,Nx,Ny);

%CALCULATE SIGNAL AND PHASE FOR EACH PIXEL AND CORRESPONDING FREQUENCY
for s=1000:ns,  %exclude dummy scans
for i=1:Nx,
    for j=1:Ny,
        if ph1(i,j)~=0
        freq=ceil(freqoff(i,j))+501;
        sig = sqrt((signal(freq,1,s)+(signal(freq,2,s))).^2);
        si = signal(freq,1,s)+sqrt(-1)*signal(freq,2,s);
        phase = angle(si);
        phantom((s-999),i,j)=ph1(i,j) * sig *(cos(phase) + sqrt(-1) * sin(phase)); %*mag(freq)*phase(freq);  % the signal and phase of ph1 w
        end
    end
    i
end

    s
end

Np = ns-1000;

%GENERATE GOLDEN ANGLE ARRAY
ang=0; 
for n =1:Np, angles(n)=ang+(n*111.246);end,

%GENERATE RADIAL PROJECTIONS USING SIGNAL AND ANGLES
a=angles;
for s=1:ns-1000
radonr(:,s)=radon(real(squeeze(phantom(s,:,:))),a(1,s)); 
radoni(:,s)=radon(imag(squeeze(phantom(s,:,:))),a(1,s));
end
imager=iradon(radonr,angles,Nx);
imagei=iradon(radoni,angles,Nx);
newim=complex(imager,imagei);
allims = newim;

%PLOT TO CHECK THE RECONSTRUCTIONS
% figure, imagesc(squeeze(abs(radonx)));
% title('fbp of data phase cycling 1');
% figure, imagesc(squeeze((imager)));
% figure, imagesc(squeeze((imagei)));
% figure, imagesc(angle((newim)));
% title('Phase of Image')

%SCALED IMAGE FROM 0 TO 1
% scaled = newim/max(max(newim));
% artefact = (1-ph1).*scaled; %APPLY INVERSE ROI MASK TO GET ARTEFACT
% figure, colormap(gray); imagesc(abs((newim)/max(max(newim))));
% figure, colormap(gray); imagesc(abs(artefact));
% inc = 360/totp;
% title(sprintf( "%s%d%s%d", 'DYPR Projections = ', ns-1000, 'Increment = ', inc))
% title(sprintf( "%s%d", 'Phase Cyling = ', p))

%SPECIFIY A SIGNAL MASK AND CALCULATE CHARACTERISTICS
% signal_filter = ph1;
% noise_filter = 1-signal_filter;
    signal_filter(signal_filter == 0) = NaN;
    shnigal = scaled.*signal_filter;
    sig_mean1 = mean(shnigal(~isnan(shnigal)));
    sig_max = max(shnigal(~isnan(shnigal)));
    sig_min = min(shnigal(~isnan(shnigal)));
    sig_total = sum(shnigal(~isnan(shnigal)));
    sig_sd = std(shnigal(~isnan(shnigal)));
%     save signal_filter signal_filter
    
%SPECIFY A NOISE/ARTEFACT MASK AND CALCULATE CHARACTERISTICS
%     noise = roipoly;
% load noise_filter
%     invnoise_filter = double(noise);
%     noise_filter = 1-invnoise_filter;
    noise_filter(noise_filter == 0) = NaN;
    naise = scaled.*noise_filter;
    noi_mean1 = mean(naise(~isnan(naise)));
    noi_max = max(naise(~isnan(naise)));
    noi_min = min(naise(~isnan(naise)));
    noi_total = sum(naise(~isnan(naise)));
    noi_sd = std(naise(~isnan(naise)));
    
%     save noise_filter noise_filter


%PRINT OUT MEASUREMENTS
words = "The max signal is: %0.5f  \nThe min signal is: %0.5f \nThe mean signal  is: %0.5f \nThe standard deviation of signal is: %0.5f \nThe total signal is: %0.5f ";
words2 ="The max artefact is: %0.5f \nThe min artefact is: %0.5f \nThe mean artefact  is: %0.5f  \nThe standard deviation of artefact is: %0.5f \nThe total artefact is: %0.5f";
sprintf(words,abs(sig_max),abs(sig_min),abs(sig_mean1),abs(sig_sd),abs(sig_total))
sprintf(words2, abs(noi_max), abs(noi_min),abs(noi_mean1), abs(noi_sd), abs(noi_total))





%% FUNCTIONS USED

%rotation matrix is for rotation about X axis
function R = rotx(ang)
ang = ang/180*pi;
R = [1, 0, 0; 0, cos(ang), sin(ang); 0, -sin(ang), cos(ang)];
end

%rotation matrix is for rotation about Z axis
function R = rotz(phi)
phi = phi/180*pi;
R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0, 0 1];
end

%rotation matrix is for rotation about Y axis
function Ry=roty(alpha)
alpha = alpha/180*pi;
Ry = [cos(alpha) 0 sin(alpha);0 1 0;-sin(alpha) 0 cos(alpha)];
end

% function [Afp,Bfp]=freeprecess(TR,bloodT1,bloodT2, myoT1, myo, df)
%  phi = 2*pi*df*T/1000;	% Resonant precession, radians.
%  E1 = exp(-TR/T1);	
%  E2 = exp(-TR/T2);
% 
%  Afp = [E2 0 0;0 E2 0;0 0 E1]*rotz(phi);
%  Bfp = [0 0 1-E1]';
% end

%HEURISTICALLY CALCULATE SSFP SIGNAL
function [magnetization]=dypr(TR,df, e1m, e2m, alpha, idx, numTR, mm, t, myoT1, myoT2, ns,f, ph)
k = 1;

numTR=ns;

    while k < 2

 
    % alpha/2 flip about x axis in clockwise direction for Mx,y,z
    % starts entering data at idx = 2
    mm(:,idx) = roty(-alpha/2) * mm(:,idx-1);
    t(idx) = t(idx-1); 
    idx = idx+1;
    
    %PRECESS AND DECAY/RELAX FOR TR/2
        phi2 = 360*df*(TR/2)/1000;
        e2m2 = exp(-(TR/2)/myoT2);       %standard exponential terms 
        e1m2 = exp(-(TR/2)/myoT1);
        Pm2 = [e2m2, 0, 0; 0, e2m2, 0; 0, 0, e1m2]*rotz(phi2);
        mm(:,idx) = Pm2 * mm(:,idx-1) + (1-e1m2) * [0;0;1];
        t(idx) = t(idx-1) + TR/2; 
        idx = idx+1; 
% %     

%ALPHA RF PULSES WITH ASSOCIATED PHASE, FOLLOWED BY TR
%PRECESSION/DECAY/RELAXATION
    for idx2 = 1:numTR
        % precession and alpha flip completed for each TR during acquisition
        % precession - rotation matrix? phi - TR or acq?
        % at index 3, repeated for the number TRs
     
        % alpha - rotate Mx,y,z by alpha about x axis
        if mod(idx2,2)==0 
            mm(:,idx) = roty(alpha) * rotz(ph(idx2)) * mm(:,idx-1);
            t(idx) = t(idx-1); 
            idx = idx+1; 
            
        else
            mm(:,idx) = roty(alpha) * rotz(ph(idx2)) * mm(:,idx-1);
            t(idx) = t(idx-1); 
            idx = idx+1; 
        end
        
        %phi = 360*df*TR/1000;
        %Pm = [e2m, 0, 0; 0, e2m, 0; 0, 0, e1m]*rotz(phi);
        mm(:,idx) = Pm2 * mm(:,idx-1) + (1-e1m2) * [0;0;1];
        t(idx) = t(idx-1) + TR/2; %time = t(3)+ 3ms
        idx = idx+1; 
        mm(:,idx) = Pm2*mm(:,idx-1) + (1-e1m2) * [0;0;1];
        t(idx) = t(idx-1) + TR/2; %time = t(3)+ 3ms
        idx = idx+1; 

        %READOUT MAGNETISATION AT TE (TR/2)
        magnetization(f, 1,idx2) = mm(1,idx-2);
        magnetization(f, 2,idx2) = mm(2, idx-2);
        
    end
   
    k = k + 1;

end
end
