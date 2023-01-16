% Recon of 4-pass radial lc-ssfp using grid toolbox
% jie xiang @yale mrrc
% can be easily adjusted to n-pass lc-ssfp
% can be followed by parallel imaging recon, eg.NLINV using BART

clear all

%SET INITIAL PARAMTERS 
res=192;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=2;    %oversampling-factor gridsize/FOV/resolution

close all;
%% IS DATA INTERLEAVED OR NOT?
isDynamicByAverage = 0; % 1 = average(acquisition), 0 = measurement(repetition)
ManualCutPhase = 20; %23 / 27
isdynamicinterleaved = 1; %1 = on, 0 = off
isphaseinterleaved = 0; % interleaved among different phases = 1
RetroUnderSam = 0; % 1 = retrospective undersampling
%if measfileIdx >23, isdynamicinterleaved=2; end

order = 1; % standard order = 1, complementary order = 2
iscrop = 1; % crop = 1, else use the original k-data
isPhaseCor = 0; % first and zero order phase correction
isPreWhitening = 0; % pre-whitening = 1
isCoilCutManual = 0; % remove bad coils by hand = 1, else use all coils
isCoilCutCROP = 0; % remove bad coils with CR+OP(outlier pruning) method
isgridfirst = 0; % grid before combine =1, else combine first
isAdapCC = 0; % acc=1, sos=0
isavifile = 0; % save = 1
issave = 0; 
issave4 = 0; % save = 1
isinterpolate = 0;

%NOTES: DONT FORGET INTERLEAVING oN OFF AND ALSO ORDER COMPLEMENTARY OR
%measfile0{1} = 'F:\CVCoding\gridding_lcssfp\MeasData\20220303_SupumLee\Chenxi\meas_MID00153_FID140712_8D_interPha_PD60.dat'; %Jie
%measfile0{1} = 'F:\CVCoding\gridding_lcssfp\MeasData\20220317_phantom_flow\meas_MID00126_FID142122_4D_BandsTest.dat'; %Jie
%measfile0{1} = 'F:\CVCoding\gridding_lcssfp\MeasData\20220409_FlowPhantomPacemaker\meas_MID00178_FID143962_4D_interleaved_PD60.dat'; %Jie
%measfile0{1} = 'F:\CVCoding\gridding_lcssfp\MeasData\20220419_HoraceZ\meas_MID00041_FID144824_4D_interleaved_PD60.dat'; %Jie
measfile0{1} = 'F:\CVCoding\gridding_lcssfp\MeasData\20211104_RongtaoJiang_metal\meas_MID00053_FID130672_jx_10phase_90_4dyn_inter.dat'; %Jie

combMethod = 'CS-SSFP';
%combMethod = 'SOS-SSFP';
%combMethod = 'MS-SSFP';
%combMethod = 'MIP-SSFP';

% for measfileIdx = 1:1%1:9
%     measfile = measfile0{measfileIdx};
% endcl
measfile = measfile0{1};
%% READ IN RAW TWIX DATA FILE
% rawdata = MultiRAID(2).rawdata;
dispopt = 'on';
%     for measure = 1:MultiRAID(2).sMDH.sLC.ushAcquisition 
%         newraw=zeros(10,MultiRAID(2).sMDH.sLC.ulChannelId,MultiRAID(2).sMDH.sLC.ushPhase,60,MultiRAID(2).sMDH.ushSamplesInScan); %%2D
%         for phase=1:MultiRAID(2).sMDH.sLC.ushPhase %10
%             for coil=1:MultiRAID(2).sMDH.sLC.ulChannelId %26
%                 for lc = 1:60
%                     for slice = 1:10
%                     lin_idx = find(MultiRAID(2).loopcounters(:,2) == slice & MultiRAID(2).loopcounters(:,1) == lc & MultiRAID(2).loopcounters(:,2) == measure & MultiRAID(2).loopcounters(:,6) == phase & MultiRAID(2).loopcounters(:,15) == coil & MultiRAID(2).loopcounters(:,16) == MultiRAID(2).sMDH.ushSamplesInScan);                     
%                     if numel(lin_idx)>=1
%                         lin_idx = lin_idx(1);
%                         newraw(slice,coil,phase,lc,:)=rawdata(lin_idx,:);
%                     end
%                     end
%                 end 
%             end    
%         end
% %%!OPTIONAL IF READING OUT MORE THAN ONE MATRIX! jsr
%         if measure == 1; first = newraw;end
%         if measure == 2; second = newraw;end
%         if measure == 3; third = newraw;end
%         if measure == 4; fourth = newraw;end
%     end
% Sorting DATA
[ first, second, third, fourth, rawdata, ~ ,MultiRAID, ~,~ ] = JX_ReadSiemensMeasVD13_2DLCSSFP(measfile,isDynamicByAverage, dispopt) ;
% [NoCoil,NoPhase,Np,Nx] = size(first);
%third = first*0;
%fourth = first*0;
%%
preWflag = 0;
if size(MultiRAID,2)==2
    phases = MultiRAID(2).sMDH.sLC.ushPhase; %PHASES
    coils= MultiRAID(2).sMDH.ulChannelId;%COILS
    if isDynamicByAverage == 1
        measure= MultiRAID(2).sMDH.sLC.ushAcquisition;
    elseif isDynamicByAverage == 0
        measure= MultiRAID(2).sMDH.sLC.ushRepetition; %MEASURES / PASSES
    end
    part = MultiRAID(2).sMDH.sLC.ushPartition;   %SLICES
    Np = MultiRAID(2).sMDH.sLC.ushLine;      %PROJECTIONS PER PASS
    Nx = MultiRAID(2).sMDH.ushSamplesInScan; %SAMPLES
    if isPreWhitening == 1
        tic;
        disp('----------------- pre Whitening -------------------');
        [first,rn] = LcssfpPreWhiten(MultiRAID(1).rawdata,coils,first);
        [second,~] = LcssfpPreWhiten(MultiRAID(1).rawdata,coils,second);
        [third,~] = LcssfpPreWhiten(MultiRAID(1).rawdata,coils,third);
        [fourth,~] = LcssfpPreWhiten(MultiRAID(1).rawdata,coils,fourth);
        preWflag = 1;
        toc;
    end
end

if size(MultiRAID,2)==1
    phases = MultiRAID.sMDH.sLC.ushPhase; %PHASES
    coils= MultiRAID.sMDH.ulChannelId;%COILS
    if isDynamicByAverage == 1
        measure= MultiRAID.sMDH.sLC.ushAcquisition;
    elseif isDynamicByAverage == 0
        measure= MultiRAID.sMDH.sLC.ushRepetition; %MEASURES / PASSES
    end
    part = MultiRAID.sMDH.sLC.ushPartition;   %SLICES
    Np = MultiRAID.sMDH.sLC.ushLine;      %PROJECTIONS PER PASS
    Nx = MultiRAID.sMDH.ushSamplesInScan; %SAMPLES
end
%clear MultiRAID
if isDynamicByAverage == 1
    %Manually choose
    %phases = ManualCutPhase;
    if isinterpolate == 1
    disp('----------------- retro gating interpolating -------------------');
    figure(101)
    subplot 211
    imagesc(abs(squeeze(first(1,:,:,1))))
    parfor k  = 1:Nx
        raw(:,:,:,k) = retro_gating_fill(first(:,:,:,k), 2) ;
    end
    subplot 212
    imagesc(abs(squeeze(raw(1,:,:,1))))
    first = raw ;
    
    parfor k  = 1:Nx
        raw(:,:,:,k) = retro_gating_fill(second(:,:,:,k), 2) ;
    end
    second = raw ;
    
    parfor k  = 1:Nx
        raw(:,:,:,k) = retro_gating_fill(third(:,:,:,k), 2) ;
    end
    third = raw ;
    
    parfor k  = 1:Nx
        raw(:,:,:,k) = retro_gating_fill(fourth(:,:,:,k), 2) ;
    end
    fourth = raw ;
    
    else
    %Manually choose
    phases = ManualCutPhase;
    first = first(:,1:phases,:,:);
    second = second(:,1:phases,:,:);
    third = third(:,1:phases,:,:);
    fourth = fourth(:,1:phases,:,:);
    end
end

%GENERATE PROJECTION ANGLES DEPENDING ON ODD OR EVEN PROJECTIONS
if mod(Np,2) == 0
    phival=[0:pi/Np:pi-pi/Np];
end
if mod(Np,2) ~=0
    phival=[0:2*pi/Np:2*pi-2*pi/Np];
end
%TAKE OUTPUT OF READIN FUNCTION AND ORGANISE INTO DESIRED FORMAT
nunder=1;% undersampling factor
% phival= phival(1:nunder:end);
if RetroUnderSam == 1
   tic;
   disp('----------------- retrospective undersampling -------------------');
   Np=Np/4;
   first = LcssfpUnderSam(first);
   second = LcssfpUnderSam(second);
   third = LcssfpUnderSam(third);
   fourth = LcssfpUnderSam(fourth);
   phival=[0:pi/Np:pi-pi/Np];
   toc;
end
ksampsone = first;
ksampsnewone=permute(ksampsone,[1 2 4 3]);
ksampsnewone= ksampsnewone(:,:,:,:,1:nunder:end);
ksampstwo = second;
ksampsnewtwo=permute(ksampstwo,[1 2 4 3]);
ksampsnewtwo= ksampsnewtwo(:,:,:,:,1:nunder:end);
ksampsthree = third;
ksampsnewthree=permute(ksampsthree,[1 2 4 3]);
ksampsnewthree= ksampsnewthree(:,:,:,:,1:nunder:end);
ksampsfour = fourth;
ksampsnewfour=permute(ksampsfour,[1 2 4 3]);
ksampsnewfour= ksampsnewfour(:,:,:,:,1:nunder:end);


%% JX coil selection(manual)
% bad coils 5,7,8.20,21,29
if isCoilCutManual == 1
disp('----------------- get rid of the bad coils -------------------');
tic;
ksampsnewone = ksampsnewone([2,4,5,6,7,9,10,13,14,18,25,26,27,28,31,33,34],:,:,:);
ksampsnewtwo = ksampsnewtwo([2,4,5,6,7,9,10,13,14,18,25,26,27,28,31,33,34],:,:,:);
ksampsnewthree = ksampsnewthree([2,4,5,6,7,9,10,13,14,18,25,26,27,28,31,33,34],:,:,:);
ksampsnewfour = ksampsnewfour([2,4,5,6,7,9,10,13,14,18,25,26,27,28,31,33,34],:,:,:);
coils=size(ksampsnewone,1);
toc;
end
   

%% JX crop in image domain
if iscrop ==1
disp('------------------- croping the data ------------------------');
tic;
Ksampsnewone = zeros(coils,phases,Nx/oversmpl,Np);
for nphases = 1:phases
    for ncoils = 1:coils
        Ktemp1 = fftshift(fft(fftshift(squeeze(ksampsnewone(ncoils,nphases,:,:))),[],1));
        Ktemp2 = Ktemp1(Nx/4+1:Nx*3/4,:);
        Ktemp3 = fftshift(ifft(fftshift(Ktemp2)));
        Ksampsnewone(ncoils,nphases,:,:) = Ktemp3;
    end
end

Ksampsnewtwo = zeros(coils,phases,Nx/oversmpl,Np);
for nphases = 1:phases
    for ncoils = 1:coils
        Ktemp1 = fftshift(fft(fftshift(squeeze(ksampsnewtwo(ncoils,nphases,:,:))),[],1));
        Ktemp2 = Ktemp1(Nx/4+1:Nx*3/4,:);
        Ktemp3 = fftshift(ifft(fftshift(Ktemp2)));
        Ksampsnewtwo(ncoils,nphases,:,:) = Ktemp3;
    end
end

Ksampsnewthree = zeros(coils,phases,Nx/oversmpl,Np);
for nphases = 1:phases
    for ncoils = 1:coils
        Ktemp1 = fftshift(fft(fftshift(squeeze(ksampsnewthree(ncoils,nphases,:,:))),[],1));
        Ktemp2 = Ktemp1(Nx/4+1:Nx*3/4,:);
        Ktemp3 = fftshift(ifft(fftshift(Ktemp2)));
        Ksampsnewthree(ncoils,nphases,:,:) = Ktemp3;
    end
end

Ksampsnewfour = zeros(coils,phases,Nx/oversmpl,Np);
for nphases = 1:phases
    for ncoils = 1:coils
        Ktemp1 = fftshift(fft(fftshift(squeeze(ksampsnewfour(ncoils,nphases,:,:))),[],1));
        Ktemp2 = Ktemp1(Nx/4+1:Nx*3/4,:);
        Ktemp3 = fftshift(ifft(fftshift(Ktemp2)));
        Ksampsnewfour(ncoils,nphases,:,:) = Ktemp3;
    end
end
ksampsnewone = Ksampsnewone;
ksampsnewtwo = Ksampsnewtwo;
ksampsnewthree = Ksampsnewthree;
ksampsnewfour = Ksampsnewfour;
clear Ksampsnewone Ksampsnewtwo Ksampsnewthree Ksampsnewfour Ktemp1 Ktemp2 Ktemp3
toc;
Nx = Nx/oversmpl;
oversmpl = 1;
end



%% JX first and zero order phase correction
if isPhaseCor ==1
disp('----------------- Phase Correction ------------------');
tic;
%Ksampsnewone = ksampsnewone;
for nphases = 1:phases
    for ncoils = 1:coils
        ToCor = squeeze(ksampsnewone(ncoils,nphases,:,:));
        ksampsnewone(ncoils,nphases,:,:) = LcssfpPhaseCor(ToCor);
    end
end

for nphases = 1:phases
    for ncoils = 1:coils
        ToCor = squeeze(ksampsnewtwo(ncoils,nphases,:,:));
        ksampsnewtwo(ncoils,nphases,:,:) = LcssfpPhaseCor(ToCor);
    end
end

for nphases = 1:phases
    for ncoils = 1:coils
        ToCor = squeeze(ksampsnewthree(ncoils,nphases,:,:));
        ksampsnewthree(ncoils,nphases,:,:) = LcssfpPhaseCor(ToCor);
    end
end

for nphases = 1:phases
    for ncoils = 1:coils
        ToCor = squeeze(ksampsnewfour(ncoils,nphases,:,:));
        ksampsnewfour(ncoils,nphases,:,:) = LcssfpPhaseCor(ToCor);
    end
end
clear ToCor
toc;
end



%% reconstruction
disp('------------------- Reconstruction starts -------------------------');
% density compensation filter...
% create Shepp-Logan dcf (looks smoother than ramlak)
tic;
disp('---------- data preparation and linear combination ------------');
tic;
if isgridfirst == 1
    k=[-.5:1/Nx:.5-1/Nx];
    dcf1=abs(k).*cos(k);
    dcf=repmat(dcf1,[Np 1]); 

    %INITIALIZE TO SAVE COMPUTATION TIME
    ft_LC_0 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    ft_LC_90 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    ft_LC_180 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    ft_LC_270 = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    add = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    sub = zeros(phases, coils, oversmpl*res,oversmpl*res); 
    four = zeros(phases, coils, oversmpl*res,oversmpl*res); 

    %create dcf-matrix (Ram-lak-filter) (same for all coildata)
    % dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
    % % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
    % dcf2=dcf1/max(dcf1);
    % dcf=repmat(dcf2,[Np 1]); %384*39 matrix
    % figure;plot(dcf(1,:));
    % end % if correct_all
    % locate kspace samples in complex space (same for all coildata)

    %GENERATE TRAJECTORIES FOR GRIDDING
    loc=zeros(measure,phases,Nx,Np);
    if isdynamicinterleaved == 1
        for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
            for j=1:Nx 
                for k=1:Np
                    for nphase = 1:phases
                        if m == 1
                            pha=phival(k);  %NO OFFSET FOR FIRST PASS
                        end
                        if m > 1 && mod(Np,2) ==0
                            offset = (pi)/(measure*Np);    %SPECIFY OFFSET FOR EVEN Np
                        end
                        if m >1 && mod(Np,2) ~=0
                            offset = (2*pi)/(measure*Np);  %SPECIFY OFFSET FOR ODD Np
                        end
                        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
                        if m == 2
                            pha=phival(k)+(2*offset);
                        end
                        if m == 3
                            pha=phival(k)+(1*offset);
                        end
                        if m == 4
                            pha=phival(k)+(3*offset);
                        end
                    
                        %pha=phival(k)+((m-1)*offset);
                        %loc(m,j,k)=mag*exp(1i*pha);
                    
                        if  isphaseinterleaved == 1 && mod(nphase,2) == 0
                            loc(m,nphase,j,k)=mag*exp(1i*(pha+(pi)/(2*Np)));
                        else
                            loc(m,nphase,j,k)=mag*exp(1i*pha);
                        end
                    end
                end
            end
        end
    end

    if isdynamicinterleaved==0
        for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
            for j=1:Nx
                for k=1:Np
                    for nphase = 1:phases
                        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
                        pha=phival(k);
                        %loc(m,nphase,j,k)=mag*exp(1i*pha);
                        if  isphaseinterleaved == 1 && mod(nphase,2) == 0
                            loc(m,nphase,j,k)=mag*exp(1i*(pha+(pi)/(2*Np)));
                        else
                            loc(m,nphase,j,k)=mag*exp(1i*pha);
                        end
                    end
                end 
            end
        end
    end

% FOR EACH SLICE, EACH COIL, EACH PHASE, GRID RADIAL DATA, PERFORM 2D FFT
warning('off')
for cl = 1:coils
    for ph = 1:phases
%         size(squeeze(ksampsnewone(cl,ph,:,:)))
        datone(ph,cl,:,:) = gridkb(loc(1,ph,:,:),squeeze((ksampsnewone(cl,ph,:,:))),dcf',oversmpl*res,kwidth,oversmpl);
        imone(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datone(ph,cl,:,:)))))));
%         imone(ph,cl,:,:) = deapod(ifftshift(fft2(fftshift(squeeze(datone(ph,cl,:,:))))),oversmpl,res,kwidth);
        dattwo(ph,cl,:,:) = gridkb(loc(2,ph,:,:),squeeze((ksampsnewtwo(cl,ph,:,:))),dcf',oversmpl*res,kwidth,oversmpl);
        imtwo(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(dattwo(ph,cl,:,:)))))));
        datthree(ph,cl,:,:) = gridkb(loc(3,ph,:,:),squeeze((ksampsnewthree(cl,ph,:,:))),dcf',oversmpl*res,kwidth,oversmpl);
        imthree(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datthree(ph,cl,:,:)))))));
        datfour(ph,cl,:,:) = gridkb(loc(4,ph,:,:),squeeze((ksampsnewfour(cl,ph,:,:))),dcf',oversmpl*res,kwidth,oversmpl);
        imfour(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datfour(ph,cl,:,:)))))));

    end
end 
warning('on')
% end
%
%PERFORM LINEAR COMBINATION (CAN ALSO BE DONE IN FOURIER SPACE)
%NOTE WHICH IS 90 AND WHICH IS 180!! DEPENDS ON ACQUISITION ORDER
%FOR STANDARD ORDER FT_LC_90 = IMTWO, FT_LC_180 = IMTHREE
%FOR COMPLEMENTARY ORDER FT_LC_90 = IMTHREE, FT_LC_180 = IMTWO
if order == 1
for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_90(phase, coil,:,:) = (imtwo(phase, coil, :, :)); 
        ft_LC_180(phase, coil,:,:) = (imthree(phase, coil, :, :)); 
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        if strcmp(combMethod, 'CS-SSFP')
            four(phase, coil,:,:) = abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:));
            %four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
        elseif strcmp(combMethod, 'SOS-SSFP') %SOS-SSFP
            four(phase, coil,:,:) = sqrt(abs(ft_LC_0(phase, coil,:,:)).^2 + abs(ft_LC_90(phase, coil, :,:)).^2 + abs(ft_LC_180( phase, coil,:,:)).^2 + abs(ft_LC_270( phase, coil, :,:)).^2);
        elseif strcmp(combMethod, 'MS-SSFP') %MS-SSFP
            four(phase, coil,:,:) = (abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:)));
        elseif strcmp(combMethod, 'MIP-SSFP') %MIP-SSFP
            temp(:,:,1)=squeeze(ft_LC_0(phase, coil,:,:));temp(:,:,2)=squeeze(ft_LC_90(phase, coil,:,:));temp(:,:,3)=squeeze(ft_LC_180(phase, coil,:,:));temp(:,:,4)=squeeze(ft_LC_270(phase, coil,:,:));
            four(phase, coil,:,:) = max(temp,[],3);
        end
    end
end
end
if order == 2
    for coil = 1:coils
    for phase = 1:phases
        ft_LC_0(phase, coil,:,:) = (imone(phase, coil, :, :));
        ft_LC_180(phase, coil,:,:) = (imtwo(phase, coil, :, :)); 
        ft_LC_90(phase, coil,:,:) = (imthree(phase, coil, :, :));
        ft_LC_270(phase, coil,:,:) = (imfour(phase, coil, :, :));
        add(phase, coil,:,:) = ft_LC_0( phase, coil,:,:)+((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         add(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)+((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        sub(phase, coil,:,:) = ft_LC_0(phase, coil,:,:)-((sqrt(-1))*(ft_LC_180(phase, coil, :,:)));
%         sub(phase, coil,:,:) = ft_LC_180(phase, coil,:,:)-((sqrt(-1))*(ft_LC_270(phase, coil, :,:)));
        if strcmp(combMethod, 'CS-SSFP')
            four(phase, coil,:,:) = abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:));
            %four(phase, coil,:,:) = (ft_LC_0(phase, coil,:,:)) + (ft_LC_90(phase, coil, :,:)) + (ft_LC_180( phase, coil,:,:)) + (ft_LC_270( phase, coil, :,:));
        elseif strcmp(combMethod, 'SOS-SSFP') %SOS-SSFP
            four(phase, coil,:,:) = sqrt(abs(ft_LC_0(phase, coil,:,:)).^2 + abs(ft_LC_90(phase, coil, :,:)).^2 + abs(ft_LC_180( phase, coil,:,:)).^2 + abs(ft_LC_270( phase, coil, :,:)).^2);
        elseif strcmp(combMethod, 'MS-SSFP') %MS-SSFP
            four(phase, coil,:,:) = (abs(ft_LC_0(phase, coil,:,:)) + abs(ft_LC_90(phase, coil, :,:)) + abs(ft_LC_180( phase, coil,:,:)) + abs(ft_LC_270( phase, coil, :,:)));
        elseif strcmp(combMethod, 'MIP-SSFP') %MIP-SSFP
            temp(:,:,1)=squeeze(ft_LC_0(phase, coil,:,:));temp(:,:,2)=squeeze(ft_LC_90(phase, coil,:,:));temp(:,:,3)=squeeze(ft_LC_180(phase, coil,:,:));temp(:,:,4)=squeeze(ft_LC_270(phase, coil,:,:));
            four(phase, coil,:,:) = max(temp,[],3);
        end
    end
    end
end

else % combine before the griding
    k=[-.5:1/Nx:.5-1/Nx];
    dcf1=abs(k).*cos(k);
    dcf=repmat(dcf1,[measure*Np 1]); 

    loc=zeros(measure,phases,Nx,Np);  
    location =zeros(Nx,phases,measure*Np);

    if isdynamicinterleaved == 1
    for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
    for j=1:Nx 
        for k=1:Np
            for nphase = 1:phases
                if m == 1
                     pha=phival(k);  %NO OFFSET FOR FIRST PASS
                end
                if m > 1 && mod(Np,2) ==0
                    offset = (pi)/(measure*Np);    %SPECIFY OFFSET FOR EVEN Np
                end
                if m >1 && mod(Np,2) ~=0
                    offset = (2*pi)/(measure*Np);  %SPECIFY OFFSET FOR ODD Np
                end
                mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
                if m == 2
                    pha=phival(k)+(2*offset);
                end
                if m == 3
                    pha=phival(k)+(1*offset);
                end
                if m == 4
                    pha=phival(k)+(3*offset);
                end
            %pha=phival(k)+((m-1)*offset);
            %loc(m,j,k)=mag*exp(1i*pha);
                if  isphaseinterleaved == 1 && mod(nphase,2) == 0
                    loc(m,nphase,j,k)=mag*exp(1i*(pha+(pi)/(2*Np)));
                else
                    loc(m,nphase,j,k)=mag*exp(1i*pha);
                end
            end
        end
    end
    end
    end

    if isdynamicinterleaved==0
        for m =1:measure  %LOOP THROUGH PASSES / CALCULATE EACH SET OF TRAJECTORIES
            for j=1:Nx
                for k=1:Np
                    for nphase = 1:phases
                        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
                        pha=phival(k);
                        %loc(m,j,k)=mag*exp(1i*pha);
                        if  isphaseinterleaved == 1 && mod(nphase,2) == 0
                            loc(m,nphase,j,k)=mag*exp(1i*(pha+(pi)/(2*Np)));
                        else
                            loc(m,nphase,j,k)=mag*exp(1i*pha);
                        end
                    end
                end 
            end
        end
    end

    for lpmea = 1:measure
        for nphase = 1:phases
            location(:,nphase,((lpmea-1)*Np+1):lpmea*Np)=loc(lpmea,nphase,:,:);
        end
    end

    for cl = 1:coils
        for ph = 1:phases
            for m = 1:measure
                if  m == 1
                    ksampsnew(cl,ph,:,((m-1)*Np+1):m*Np)=ksampsnewone(cl,ph,:,:);
                elseif m == 2
                    ksampsnew(cl,ph,:,((m-1)*Np+1):m*Np)=ksampsnewtwo(cl,ph,:,:);
                elseif m == 3
                    ksampsnew(cl,ph,:,((m-1)*Np+1):m*Np)=ksampsnewthree(cl,ph,:,:);
                else    
                    ksampsnew(cl,ph,:,((m-1)*Np+1):m*Np)=ksampsnewfour(cl,ph,:,:);
                end
            end
            datfour(ph,cl,:,:) = gridkb(location(:,ph,:,:),squeeze((ksampsnew(cl,ph,:,:))),dcf',oversmpl*res,kwidth,oversmpl);
            four(ph,cl,:,:) = (ifftshift(fft2(fftshift((squeeze(datfour(ph,cl,:,:)))))));
        end
    end 
end
toc;

if isCoilCutCROP == 1
disp('-------- get rid of the bad coils with CR + OP -----------');
tic;

FilterBlock = zeros(Nx,Nx);
FilterBlock((Nx/2-30):(Nx/2+30),(Nx/2-30):(Nx/2+30)) = 1;
for ncoil = 1:coils
    for nphase = 1:phases
        Ij = squeeze(four(nphase,ncoil,:,:));
        IjTemp = ifftshift(ifft2(fftshift(Ij))).*FilterBlock;
        Lj = fftshift(fft2(ifftshift(IjTemp)));
        DeltaIL = Ij - Lj;
        SR(ncoil,nphase) = norm(DeltaIL,2)/norm(Lj,2); % ||Ij-Lj||2 / ||Lj||2
    end
end
SR = sqrt(sum(SR.^2,2))./phases;
CoilIdx = kmeans(SR,2);
SRthreshold = 0.2;

cflag = 1;
SRthreshold = 0.15;
for ncoil = 1:coils
    if SR(ncoil) < SRthreshold
        fourtemp(:,cflag,:,:) = four(:,ncoil,:,:);
        fourtemp2(:,cflag,:,:) = datfour(:,ncoil,:,:);
        %goodcoils(cflag) = ncoil;
        cflag = cflag + 1;
    end
end
four = fourtemp;
datfour = fourtemp2;
coils = cflag - 1;
clear fourtemp


if mean(SR(CoilIdx==2))- mean(SR(CoilIdx==1)) > 2*std(SR)
    cflag = 1;
    for ncoil = 1:coils
        % calculate the SR value(streaking artifact ratio)
        if CoilIdx(ncoil) == 1
            fourtemp(:,cflag,:,:) = four(:,ncoil,:,:);
            fourtemp2(:,cflag,:,:) = datfour(:,ncoil,:,:);
            %goodcoils(cflag) = ncoil;
            cflag = cflag + 1;
        end
    end
    four = fourtemp;
    datfour = fourtemp2;
    coils = cflag - 1;
    clear fourtemp
end
toc;
end 


if isphaseinterleaved == 1
tic;
disp('----------- interleaved phase UNFOLD ------------');
% CropWindow = hann(phases); % hanning window, too strong
n=32;
t=-n:n;
filter=1./(1+exp(t-30/1))-1./(1+exp(t+30/1));
CropWindow=resample(filter,double(phases-1),64); %fermi filter
% plot(CropWindow);
for ncoil = 1:coils
    for xi = 1:Nx
        for yi = 1:Nx
            TimeCurve(:,ncoil,xi,yi) = four(:,ncoil,xi,yi);
            % TimeCurve = abs(four(:,ncoil,xi,yi));
            FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,ncoil,xi,yi))))).*(CropWindow');
            four(:,ncoil,xi,yi)=fftshift(ifft(ifftshift(FreCurve)));
            %four(:,ncoil,xi,yi)=FreCurve;
        end
    end
end
toc;
end

if isAdapCC == 1
disp('----------- combine different coils with ACC ------------');
tic;
four2=zeros(Nx,Nx,phases);
parfor ph = 1:phases
    if preWflag == 1
        four2(:,:,ph) = AdaptiveCombine(squeeze(four(ph,:,:,:)),1,rn);
    else
        four2(:,:,ph) = AdaptiveCombine(squeeze(four(ph,:,:,:)),1);
    end
end
toc;
else
disp('----------- combine different coils with SOS ------------');
tic;
four2 = squeeze(sqrt(sum(abs(four).^2,2)));
four2 = permute(four2,[2,3,1]);
toc;
end
close all;

brightening = 0.3;
windowing   = 0.8;
figure; 
for idx = 1:phases
img_disp = imresize(squeeze(abs(four2(:,:,idx))),2);
imshow(brighten(img_disp/max(img_disp(:)),brightening), [0 windowing]);title(['phase = ',num2str(idx)]), pause(0.2)
end

figure,
for idx = 1:phases
img_disp = imresize(squeeze(abs(four2(:,:,idx))),2);
imshow(img_disp, [0,0.25]); pause(0.2)
end

%%
figure,idx = 9;
windowing   = 0.9; brightening = 0.01;
img_disp = imresize(squeeze(abs(four2(49:144,62:142,idx))),1.5);
imshow((fliplr(rot90(brighten(img_disp/max(img_disp(:)),brightening),3))), [0 windowing]);title(['phase = ',num2str(idx)])
%%

figure,imshow(imresize(squeeze(abs(four2(:,:,idx))),2),[])
%% save
if issave == 1
disp('------------------- Save the images and gif -------------------------');
tic;
% save('Recon_Image\Dana\DP39_4d_pd60_20p_acc_ave.mat','four2','datfour');
save('Recon_Image\rongtao1104\RT53_0621_GRID_CR.mat','four2','datfour','four');

toc;
end

if issave4 == 1
disp('-------------- Save the 4 separate images and gif -----------------');
tic;
LcssfpGen4Img;
toc;
end


