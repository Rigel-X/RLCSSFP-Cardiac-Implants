% 2022/10/5, simulate lc-ssfp with and without interleaving
% use dicom image as base image, instead of Shepp Logan phantom, as per
% reviewer2 @JMRI, manuscript revision submitted 2022/10/31, accepted 11/1
% radial kspace created using @BART, superimposed bssfp banding mask

% different gaussian noise level and projections are simulated
%     -----  GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
%     -----  Nplist = [24,48,60,80,96,120,140,160,192];

% jie xiang @yale mrrc

clear all
close all
clc

%% read cine image data
addpath('F:\CVData\phase contrast\20220706_QL\MRRC_RESEARCH_JEROME_20220706_173436_244000\TRUFI_SAX_CHAMBER_CINE_0017') 
files=dir('F:\CVData\phase contrast\20220706_QL\MRRC_RESEARCH_JEROME_20220706_173436_244000\TRUFI_SAX_CHAMBER_CINE_0017\*.IMA');

rawdata = zeros(192,192,length(files));
for j=1:length(files)
rawdata(:,:,j)=double(dicomread(files(j).name));
end
rawdata = rawdata(30:144,30:144,7:16:end);

CineImg = zeros(192,192,30);
figure,
for j = 1:30
    CineImg(:,:,j) = imresize(squeeze(rawdata(:,:,j)),192/115);
    imshow(brighten(abs(squeeze(CineImg(:,:,j))),0.4),[])
end
[Nx,Ny,NoPhase] = size(CineImg); 
clear rawdata files j
save('lcssfp_simudata.mat','CineImg')

%% bssfp response profile
T1=1800;T2=200;%water
TR=3.5; % ms
freqlim=2100/TR; %Hz
alpha=35*pi/180; % optimal flip angle = arccos(exp(-TR/T1))

% different phase offsets
phaseoffset = [0,0.25,0.5,0.75]; % four dynamics
% phaseoffset = [0,1/3,2/3]; % three dynamics
% phaseoffset = [0,1/2]; % two dynamics
% phaseoffset = [0]; % one dynamic
StepHz = 0.01;

E1=exp(-TR/T1); E2=exp(-TR/T2);
a=-(1-E1)*E2*sin(alpha);  b=(1-E1)*sin(alpha); c=E2*(E1-1)*(1+cos(alpha)); d=1-E1*cos(alpha)-(E1-cos(alpha))*E2*E2;
freq1=(1/(TR*.001)).*phaseoffset;% shift in frequency, to demonstrate phase offset. 
j=1;
phasecycleNumber = 1;
freqi = zeros(1,1+2*freqlim);
for freq=-freqlim:StepHz:freqlim %hz
    phi=2*pi*(freq+freq1(phasecycleNumber))*TR*.001;
    Mssxy(j,:)=abs((a.*exp(-1i*phi)+b)./(c*cos(phi)+d));
    freqi(j)=freq;
    j=j+1;
end
figure,plot(freqi,abs(Mssxy),'r-');

%% generate intensity mask
% dephasing mask
PmkX = 32; PmkY = 82; % position of the pacemaker
DephMask = ones(192,192);
for i = 1:192
    for j = 1:192
        dist = ((i-(PmkX-5))^2+(j-(PmkY-5))^2);
        if dist <= 300
            DephMask(i,j) = 0;
        elseif dist <= 400
            DephMask(i,j) = (dist-300)/100;
        end
    end
end

% banding mask
Nx = 192; Ny = 192;NoPhase = 30;
BandMask = zeros(4,Nx,Ny);
Fmap = zeros(Nx,Ny);
CinePmk = zeros(Nx,Ny,NoPhase,4);
for i = 1:Nx
    for j = 1:Ny
        for Npc = 1:4           
            dist = ((i-PmkX)^2+(j-PmkY)^2);
            Fmap(i,j) = 600*dist^(1/7);
            phi=2*pi*(Fmap(i,j)+freq1(Npc))*TR*.001;
            BandMask(Npc,i,j)=(a.*exp(-1i*phi)+b)./(c*cos(phi)+d);
            for nphases = 1:NoPhase
                CinePmk(i,j,nphases,Npc) = abs(CineImg(i,j,nphases)).*abs(BandMask(Npc,i,j)).*exp(1i*angle(BandMask(Npc,i,j))).*DephMask(i,j);
            end
        end
        test(i,j) = max(max(abs(BandMask(:,i,j))));
    end
end
figure,imshow(Fmap,[])
min(abs(BandMask(:)))
max(abs(BandMask(:)))
figure,for i = 1:4; subplot(2,2,i),imshow(squeeze(BandMask(i,:,:)),[0,max(abs(BandMask(:)))]);end
test2 = BandMask(1,:,:)+BandMask(2,:,:)+BandMask(3,:,:)+BandMask(4,:,:);
figure,subplot(1,2,1),imshow(squeeze(test),[0,4*max(abs(BandMask(:)))]),title('MIP');subplot(1,2,2),imshow(squeeze(abs(test2)),[0,4*max(abs(BandMask(:)))]),title('lc');

figure,set(gcf,'Color',[1 1 1]);
for j = 1:30
    subplot(1,2,1),imshow(brighten(abs(squeeze(CinePmk(:,:,j,1))),0.4),[]),title(['one phase cycle, phase = ',num2str(j)])
    subplot(1,2,2),imshow(brighten(abs(squeeze(CinePmk(:,:,j,1))+squeeze(CinePmk(:,:,j,2))+squeeze(CinePmk(:,:,j,3))+squeeze(CinePmk(:,:,j,4))),0.4),[]),title('lcSSFP'),pause(0.05)
    creategif = 0;
    if creategif
        nn=getframe(gcf);
        im=frame2im(nn);
        [I,map]=rgb2ind(im,256);
        if j==1
            imwrite(I,map,'simu.gif','gif','loopcount',inf,'Delaytime',0.1)
        else
            imwrite(I,map,'simu.gif','gif','writemode','append','Delaytime',0.1)
        end
    end
end
save('lcssfp_simudata.mat','CinePmk','-append')

%% generate radial kspace
% Np = 240;
% rawkfft = zeros(384,240,30,4);
% theta = (0:Np/(Np-1): Np)*360/Np;% # of angles for simulation...over 360 deg
% for nphases = 1:30
%     for Npc = 1:4
%         [sino,xp] = radon(squeeze(CinePmk(:,:,nphases,Npc)),theta,384);%this is the sinogram data.
%         rawkfft(:,:,nphases,Npc)=fftshift(fft(double(sino),[],1),1);
%     end
% end
% save('lcssfp_simudata.mat','rawkfft','-append')
% 
% first=permute(squeeze(rawkfft(:,:,:,1)),[4,3,2,1]);
% second=permute(squeeze(rawkfft(:,:,:,2)),[4,3,2,1]);
% third=permute(squeeze(rawkfft(:,:,:,3)),[4,3,2,1]);
% fourth=permute(squeeze(rawkfft(:,:,:,4)),[4,3,2,1]);

%% generate noise-free radial ksp using BART
Nplist = [24,48,60,80,96,120,140,160,192];
load('lcssfp_simudata.mat', 'CinePmk')

for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
tmp_traj = bart(['traj -x 384 -y ', num2str(4*Np),' -r']); %('traj -x 384 -y 4*Np -r');
traj = bart('scale 0.5', tmp_traj);

ksp = zeros(4,384,4*Np,30);
for i = 1:4
    for np = 1:30
        temp = squeeze(CinePmk(:,:,np,i));
        ksp(i,:,:,np) = bart('nufft -t',traj, temp);
    end
end

save(['ksp_Np',num2str(Np),'.mat'],'ksp')

end

%% add white gaussian noise
GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'.mat'])
disp(['------ksp_Np',num2str(Np),'.mat, GNnoiseSNR = ',num2str(GNnoiseSNR)])

for i = 1:4
    for np = 1:30
        for nl = 1:Np
            ksp(i,:,nl,np) = awgn(squeeze(ksp(i,:,nl,np)),GNnoiseSNR,'measured');
        end
    end
end
save(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'],'ksp')
end
end

%% final kspace and images
GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
tmp_traj = bart(['traj -x 384 -y ', num2str(4*Np),' -r']); %('traj -x 384 -y 4*Np -r');
traj = bart('scale 0.5', tmp_traj);

for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'])
disp (['----- ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat -----'])

tic
% radial kspace with different interleaving
Ksplc = zeros(4,384,Np,30);trajlc = zeros(4,3,384,Np);
Kspid = zeros(4,384,Np,30);trajid = zeros(4,3,384,Np);
Kspip = zeros(4,384,Np,30);trajip = zeros(4,3,384,Np,30);

% simple lcSSFP
for Npc = 1:4
    Ksplc(Npc,:,:,:) = ksp(Npc,:,1:4:end,:);
    trajlc(Npc,:,:,:) = traj(:,:,1:4:end);
end


% interleaved dynamic idSSFP
for Npc = 1:4
    Kspid(Npc,:,:,:) = ksp(Npc,:,Npc:4:end,:);
    trajid(Npc,:,:,:) = traj(:,:,Npc:4:end);
end

% interleaved cardiac phases
for Npc = 1:4 
    Kspip(Npc,:,:,1:2:end) = ksp(Npc,:,Npc:4:end,1:2:end);
    Kspip(Npc,:,:,2:2:end) = ksp(Npc,:,(mod(Npc+1,4)+1):4:end,2:2:end); % additional offset for even phases
    for np = 1:30
        if mod(np,2) == 0
            trajip(Npc,:,:,:,np) = traj(:,:,(mod(Npc+1,4)+1):4:end);
        else
            trajip(Npc,:,:,:,np) = traj(:,:,Npc:4:end);
        end
    end
end
toc

% % recon test
% img = [];
% for i = 1:4
%     img(i,:,:) = bart('nufft -i -t',traj(:,:,1:4:end), Ksplc(i,:,:,1));
% end
% %img(i,:,:) = bart('nufft -i -t',traj, ksp(1,:,:)+ksp(2,:,:)+ksp(3,:,:)+ksp(4,:,:));
% figure,imshow(abs(squeeze(img(1,:,:))+squeeze(img(2,:,:))+squeeze(img(3,:,:))+squeeze(img(4,:,:))),[])

% NUFFT recon using BART
Imglc = zeros(4,192,192,30);
Imgid = zeros(4,192,192,30);
Imgip = zeros(4,192,192,30);

tic
parfor nd = 1:4
    for np = 1:30
        Imglc(nd,:,:,np) = bart('nufft -i -t',squeeze(trajlc(nd,:,:,:)), Ksplc(nd,:,:,np));
        Imgid(nd,:,:,np) = bart('nufft -i -t',squeeze(trajid(nd,:,:,:)), Kspid(nd,:,:,np));
        Imgip(nd,:,:,np) = bart('nufft -i -t',squeeze(trajip(nd,:,:,:,np)), Kspip(nd,:,:,np));
    end
end
toc

tic;
disp('----------- interleaved phase UNFOLD ------------');
% CropWindow = hann(phases); % hanning window, too strong
n=32;
t=-n:n;
windpwL = 30; % 30
filter=1./(1+exp(t-windpwL/1))-1./(1+exp(t+windpwL/1));
CropWindow=resample(filter,29,64); %fermi filter
for nd = 1:4
for xi = 1:192
    for yi = 1:192
        TimeCurve(:,xi,yi) = Imgip(nd,xi,yi,:);
        FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,xi,yi))))).*(CropWindow');
        Imgip(nd,xi,yi,:)=fftshift(ifft(ifftshift(FreCurve)));
    end
end
end
toc;

figure,imshow(brighten(abs(squeeze(Imglc(1,:,:,1))),0.3),[])

lcSSFP = squeeze(sum(Imglc,1));
idSSFP = squeeze(sum(Imgid,1));
ipSSFP = squeeze(sum(Imgip,1));

figure,
for np = 1:30
subplot(1,3,1),imshow(brighten(abs(squeeze(lcSSFP(:,:,np))),0.3),[]),title(['lcSSFP, phase = ',num2str(np)])
subplot(1,3,2),imshow(brighten(abs(squeeze(idSSFP(:,:,np))),0.3),[]),title(['idSSFP, phase = ',num2str(np)])
subplot(1,3,3),imshow(brighten(abs(squeeze(ipSSFP(:,:,np))),0.3),[]),title(['ipSSFP, phase = ',num2str(np)]),pause(0.1)
end

save(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'],'lcSSFP','idSSFP','ipSSFP','-append')
end
end

%% display all the images
GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];
load('lcssfp_simudata.mat', 'CineImg')
CineImg = CineImg.*DephMask;

figure(101),set(gcf, 'Position', get(0, 'Screensize'));set(gcf,'Color',[1 1 1]);
for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'])
disp (['----- ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat -----'])

for np = 1:30
subplot(2,3,1),imshow(brighten(abs(squeeze(lcSSFP(:,:,np))),0.3),[]),title(['lcSSFP, Np = ',num2str(Np)])
subplot(2,3,2),imshow(brighten(abs(squeeze(idSSFP(:,:,np))),0.3),[]),title(['idSSFP, Noise = ',num2str(GNnoiseSNR),'dB'])
subplot(2,3,3),imshow(brighten(abs(squeeze(ipSSFP(:,:,np))),0.3),[]),title(['ipSSFP, phase = ',num2str(np)])
subplot(2,3,4),imshow(abs(squeeze(lcSSFP(:,:,np)))-abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))))*mean(mean(abs(squeeze(lcSSFP(:,:,np))))),[]),title(['lcSSFP, Np = ',num2str(Np)])
subplot(2,3,5),imshow(abs(squeeze(idSSFP(:,:,np)))-abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))))*mean(mean(abs(squeeze(idSSFP(:,:,np))))),[]),title(['idSSFP, Noise = ',num2str(GNnoiseSNR),'dB'])
subplot(2,3,6),imshow(abs(squeeze(ipSSFP(:,:,np)))-abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))))*mean(mean(abs(squeeze(ipSSFP(:,:,np))))),[]),title(['ipSSFP, phase = ',num2str(np)]),pause(0.02)
end

end
end

%% Streak Artifact Ratio, ||Ij-Lj||2 / ||Lj||2
Nx = 192;
FilterBlock = zeros(Nx,Nx);
FilterBlock((Nx/2-30):(Nx/2+30),(Nx/2-30):(Nx/2+30)) = 1;

GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

tic
SRall = zeros(length(Nplist),length(GNnoiseSNRlist),3);
for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'])
disp (['----- ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat -----'])

SR = zeros(3,30);
for np = 1:30
    
    Ij = squeeze(abs(lcSSFP(:,:,np)));
    IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
    Lj = fftshift(fft2(ifftshift(IjTemp)));
    DeltaIL = Ij - Lj; SR(1,np) = norm(DeltaIL,2)/norm(Lj,2); % lcssfp
    
    Ij = squeeze(abs(idSSFP(:,:,np)));
    IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
    Lj = fftshift(fft2(ifftshift(IjTemp)));
    DeltaIL = Ij - Lj; SR(2,np) = norm(DeltaIL,2)/norm(Lj,2); % idssfp
    
    Ij = squeeze(abs(ipSSFP(:,:,np)));
    IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
    Lj = fftshift(fft2(ifftshift(IjTemp)));
    DeltaIL = Ij - Lj; SR(3,np) = norm(DeltaIL,2)/norm(Lj,2); % ipssfp
    
end
SR = sum(SR,2);
SRall(rep2,rep1,:) = SR;
% save(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'],'SR','-append')

end
end
toc
save('AllMat.mat','SRall')

figure,plot(Nplist,squeeze(SRall(:,5,1))),hold on; plot(Nplist,squeeze(SRall(:,5,2)));plot(Nplist,squeeze(SRall(:,5,3))),title('SR w/ GN50dB, different Np')
figure,plot(GNnoiseSNRlist,squeeze(SRall(3,:,1))),hold on; plot(GNnoiseSNRlist,squeeze(SRall(3,:,2)));plot(GNnoiseSNRlist,squeeze(SRall(3,:,3))),title('SR w/ Np60, different GN')

%% RMSE
GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];
load('lcssfp_simudata.mat', 'CineImg')
CineImg = CineImg.*DephMask;

tic
RMSEall = zeros(length(Nplist),length(GNnoiseSNRlist),3);
for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'])
disp (['----- ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat -----'])

RMSE = zeros(3,30);
for np = 1:30
    
    b = squeeze(abs(lcSSFP(:,:,np)))/mean(mean(abs(squeeze(lcSSFP(:,:,np)))));
    a = abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))));
    RMSE(1,np) = sqrt(mean((a(:)-b(:)).^2)); % lcssfp
    
    b = squeeze(abs(idSSFP(:,:,np)))/mean(mean(abs(squeeze(idSSFP(:,:,np)))));
    a = abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))));
    RMSE(2,np) = sqrt(mean((a(:)-b(:)).^2)); % idssfp
    
    b = squeeze(abs(ipSSFP(:,:,np)))/mean(mean(abs(squeeze(ipSSFP(:,:,np)))));
    a = abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))));
    RMSE(3,np) = sqrt(mean((a(:)-b(:)).^2)); % ipssfp
    
end
RMSE = sum(RMSE,2)
RMSEall(rep2,rep1,:) = RMSE;
% save(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'],'RMSE','-append')

end
end
toc
save('AllMat.mat','RMSEall','-append')

figure,plot(Nplist,squeeze(RMSEall(:,5,1))),hold on; plot(Nplist,squeeze(RMSEall(:,5,2)));plot(Nplist,squeeze(RMSEall(:,5,3))),title('RMSE w/ GN50dB, different Np')
figure,plot(GNnoiseSNRlist,squeeze(RMSEall(3,:,1))),hold on; plot(GNnoiseSNRlist,squeeze(RMSEall(3,:,2)));plot(GNnoiseSNRlist,squeeze(RMSEall(3,:,3))),title('RMSE w/ Np60, different GN')
%%
load('AllMat.mat', 'RMSEall')

GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

GNn = 1; Npn = 3;


RMSEall = RMSEall./192;
clf(figure(1)),figure(1),plot(Nplist,squeeze(RMSEall(:,GNn,1)),'-*','LineWidth',2),hold on; plot(Nplist,squeeze(RMSEall(:,GNn,2)),'-*','LineWidth',2);plot(Nplist,squeeze(RMSEall(:,GNn,3)),'-*','LineWidth',2),title(['RMSE w/ GN',num2str(GNnoiseSNRlist(GNn)),'dB, different Np']),xlim([20,200])
grid on;legend('lcSSFP','idSSFP','ipSSFP')

GNCut = 9;
RMSEallnew = RMSEall(:,1:GNCut,:);GNnoiseSNRlistnew = GNnoiseSNRlist(1:GNCut);
GNnoiseSNRlistnew = 1./exp(GNnoiseSNRlistnew./10);%GNnoiseSNRlistnew = GNnoiseSNRlistnew./max(GNnoiseSNRlistnew);
clf(figure(2)),figure(2),plot(GNnoiseSNRlistnew,squeeze(RMSEallnew(Npn,:,1)),'-*','LineWidth',2),hold on; plot(GNnoiseSNRlistnew,squeeze(RMSEallnew(Npn,:,2)),'-*','LineWidth',2);plot(GNnoiseSNRlistnew,squeeze(RMSEallnew(Npn,:,3)),'-*','LineWidth',2),title(['RMSE w/ Np',num2str(Nplist(Npn)),', different GN'])
grid on;legend('lcSSFP','idSSFP','ipSSFP')




%% SNR
GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];
load('lcssfp_simudata.mat', 'CineImg')
CineImg = CineImg.*DephMask;

% a = squeeze(abs(lcSSFP(:,:,np)));
% figure,imshow(abs(a), []);
% roi = images.roi.AssistedFreehand;
% draw(roi);
% BloodMask = createMask(roi);
% roi = images.roi.AssistedFreehand;
% draw(roi);
% BgdMask = createMask(roi);
% save('AllMat.mat','BgdMask','BloodMask','-append')

load('AllMat.mat', 'BgdMask', 'BloodMask')

tic
SNRall = zeros(length(Nplist),length(GNnoiseSNRlist),3);
for rep2 = 1:length(Nplist)
Np = Nplist(rep2);
for rep1 = 1:length(GNnoiseSNRlist)
GNnoiseSNR = GNnoiseSNRlist(rep1);

load(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'])
disp (['----- ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat -----'])

SNR = zeros(3,30);
for np = 1:30
    
    a = abs(squeeze(CineImg(:,:,np)))/mean(mean(abs(squeeze(CineImg(:,:,np)))));
    
    b = squeeze(abs(lcSSFP(:,:,np)))/mean(mean(abs(squeeze(lcSSFP(:,:,np)))));
    imgs = abs(b.*BloodMask); signal  = sum(imgs(:))/sum(BloodMask(:)~=0);
    imgn = (a-b).*(BgdMask);imgn = imgn((imgn~=0));noise = std(imgn(:));
    SNR(1,np) = 10*(signal/noise); % lcssfp
    
    
    b = squeeze(abs(idSSFP(:,:,np)))/mean(mean(abs(squeeze(idSSFP(:,:,np)))));
    imgs = abs(b.*BloodMask); signal  = sum(imgs(:))/sum(BloodMask(:)~=0);
    imgn = (a-b).*(BgdMask);imgn = imgn((imgn~=0));noise = std(imgn(:));
    SNR(2,np) = 10*(signal/noise); % idssfp
    
    b = squeeze(abs(ipSSFP(:,:,np)))/mean(mean(abs(squeeze(ipSSFP(:,:,np)))));
    imgs = abs(b.*BloodMask); signal  = sum(imgs(:))/sum(BloodMask(:)~=0);
    imgn = (a-b).*(BgdMask);imgn = imgn((imgn~=0));noise = std(imgn(:));
    SNR(3,np) = 10*(signal/noise); % ipssfp
    
end
SNR = sum(SNR,2)
SNRall(rep2,rep1,:) = SNR;
% save(['ksp_Np',num2str(Np),'_GN',num2str(GNnoiseSNR),'.mat'],'SNR','-append')

end
end
toc
save('AllMat.mat','SNRall','-append')


%%
load('AllMat.mat', 'SNRall')

GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

GNn = 7; Npn = 9; %1, 3, 9

SNRFactor = (sqrt(Nplist))';

SNRall = SNRall./max(SNRall(:)).*sqrt(Nplist(end));
clf(figure(1)),figure(1),plot(Nplist,squeeze(SNRall(:,GNn,1))./SNRFactor,'-*','LineWidth',2),hold on; plot(Nplist,squeeze(SNRall(:,GNn,2))./SNRFactor,'-*','LineWidth',2);plot(Nplist,squeeze(SNRall(:,GNn,3))./SNRFactor,'-*','LineWidth',2),title(['SNR w/ GN',num2str(GNnoiseSNRlist(GNn)),'dB, different Np']),xlim([20,200])
% figure,plot(GNnoiseSNRlist,squeeze(SNRall(Npn,:,1))),hold on; plot(GNnoiseSNRlist,squeeze(SNRall(Npn,:,2)));plot(GNnoiseSNRlist,squeeze(SNRall(Npn,:,3))),title(['SNR w/ Np',num2str(Nplist(Npn)),', different GN'])
grid on;legend('lcSSFP','idSSFP','ipSSFP')

GNCut = 9;
SNRallnew = SNRall(:,1:GNCut,:)./sqrt(Nplist(Npn));GNnoiseSNRlistnew = GNnoiseSNRlist(1:GNCut);
GNnoiseSNRlistnew = 1./exp(GNnoiseSNRlistnew./10);%GNnoiseSNRlistnew = GNnoiseSNRlistnew./max(GNnoiseSNRlistnew);
clf(figure(2)),figure(2),plot(GNnoiseSNRlistnew,squeeze(SNRallnew(Npn,:,1)),'-*','LineWidth',2),hold on; plot(GNnoiseSNRlistnew,squeeze(SNRallnew(Npn,:,2)),'-*','LineWidth',2);plot(GNnoiseSNRlistnew,squeeze(SNRallnew(Npn,:,3)),'-*','LineWidth',2),title(['SNR w/ Np',num2str(Nplist(Npn)),', different GN'])
grid on;legend('lcSSFP','idSSFP','ipSSFP')



%% reference images
load('lcssfp_simudata.mat', 'CineImg', 'CinePmk', 'Fmap')
figure,imshow(brighten(cat(2,abs(squeeze(CineImg(:,:,1))),5*abs(squeeze(CinePmk(:,:,1,3)))),0.5),[])
figure,imshow(Fmap,[]),colormap('jet')

GNnoiseSNRlist = [10,15,20,25,30,40,50,60,80];
Nplist = [24,48,60,80,96,120,140,160,192];

GNn = 7; Npn = 2;

load(['ksp_Np',num2str(Nplist(Npn)),'_GN',num2str(GNnoiseSNRlist(GNn)),'.mat'])
figure(3),
for nphases = 9:9
imshow(brighten(cat(2,abs(squeeze(lcSSFP(:,:,nphases))),abs(squeeze(idSSFP(:,:,nphases))),abs(squeeze(ipSSFP(:,:,nphases)))),0.5),[]),title(['phase = ',num2str(nphases)]),pause(0.1)
end

