
% folder = 'C:\Users\Dana Peters\Documents\';
% filepath = '';

%GENERATE B0 FIELD MAP FROM 2 GRE FILES WITH DIFFERENT TEs
close all

%INPUT FIRST FILE
measfile='E:\0718_19_Sim\meas_MID389_gre5te_FID74043.dat'
%measfile='E:\0718_19_SN\sn\meas_MID348_gre5te_FID74002.dat'
%measfile='E:\0718_19_ jc\meas_MID263_gre5te_FID73917.dat'
%measfile = 'E:\0718_19_Phu\thuongvol1\meas_MID227_gre5te_FID73881.dat'

%measfile = 'C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\0418_MF\meas_MID109_gre5te_FID65842.dat'

%C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1018_IO\timtriob_20181018_161230220
%measfile = 'C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1018_IO\timtriob_20181018_161230220\meas_MID236_gre5te_FID33648.dat'

% measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1411_MM\meas_MID167_gre5te_FID37875.dat';
%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1106_JA\meas_MID58_gre5teSL2_FID36505.dat'
%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1106_JA\meas_MID45_gre5teSL1_FID36492.dat' %slice 1 is best

% measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1005_YR3\meas_MID51_gre5te_correct_FID31659.dat'
%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\3110_ES\meas_MID194_gre5te_cor_FID35771.dat'
%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1203_MB\meas_MID148_gre5te_sl2_FID40143.dat'

%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1025_BR\timtriob_20181025_064338988\meas_MID294_gre5te_ax_FID35015.dat';
dispopt = 'on';
[rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17phasenew(measfile, dispopt);
phases = sMDH.sLC.ushPhase;         %phases
nc=sMDH.sLC.ulChannelId;            %coils
npe = sMDH.sLC.ushLine;             %lines
ns = sMDH.ushSamplesInScan;         %samples

%%SPECIFY DELTA TE FOR PHASE MAPPING AND CONVERT TO SECONDS
deltate = (5-3.6)/1000;             %seconds

%PERFORM 2D FFT ON FIRST IMAGE DATA
one = zeros(nc, phases, npe, ns);
for coil = 1:nc
    for phase = 1:phases
         one(coil, phase,:,:) = fftshift(fft2(fftshift(squeeze(newraw(coil, phase, :, :)))));
    end
end


%INPUT SECOND FILE 
measfile='E:\0718_19_Sim\meas_MID388_gre3_6te_FID74042.dat'
%measfile='E:\0718_19_SN\sn\meas_MID348_gre5te_FID74002.dat'
%measfile='E:\0718_19_ jc\meas_MID264_gre3_6te_FID73918.dat'
%measfile = 'E:\0718_19_Phu\thuongvol1\meas_MID228_gre3_6te_FID73882.dat'
%measfile = 'C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\0418_MF\meas_MID110_gre3_6te_FID65843.dat'
%measfile = 'C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1018_IO\timtriob_20181018_161230220\meas_MID237_gre3_6te_FID33649.dat'

% measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1411_MM\meas_MID168_gre3_6te_FID37876.dat';




%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1106_JA\meas_MID44_gre3_6teSL1_FID36491.dat' % slice 1 is best

% measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1005_YR3\meas_MID49_gre3_6te_FID31657.dat'


%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\3110_ES\meas_MID195_gre3_6te_FID35772.dat'

%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1203_MB\meas_MID147_gre3_6te_sl2_FID40142.dat'
%measfile='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1025_BR\timtriob_20181025_064338988\meas_MID295_gre3_6te_ax_FID35016.dat'

dispopt = 'on';
[rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17phasenew(measfile, dispopt);
phases = sMDH.sLC.ushPhase;         %phases
nc=sMDH.sLC.ulChannelId;            %coils
npe = sMDH.sLC.ushLine;             %lines
ns = sMDH.ushSamplesInScan;         %samples

%INITIALIZE MATRICES
two = zeros(nc, phases, npe, ns);
diff = zeros(nc, phases, npe, ns);
freq = zeros(nc, phases, npe, ns);
pha = zeros(nc, phases, npe, ns);

%PERFORM 2D FFT ON SECOND IMAGE DATA
for coil = 1:nc
    for phase = 1:phases
        two(coil, phase,:,:) = fftshift(fft2(fftshift(squeeze(newraw(coil, phase, :, :)))));
    end
end

%SUBTRACT IMAGES AND CALCULATE THE PHASE DIFFERENCE
for coil = 1:nc
    diff(coil,1,:,:) = one(coil, 1,:,:)./two(coil,1,:,:);
%     diff(coil,1,:,:) = (one(coil,1,:,:)) - (two(coil,1,:,:));   
    pha(coil,1,:,:) = angle((diff(coil,1,:,:)));%sqrt(-1)*
end
%
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1]

%CALCULATING FREQUENCY
%PLOT FREQUENCY AND PHASE MAP FOR EACH COIL
for coil = 1:nc
    freq(coil,phase,:,:) = ((pha(coil,phase,:,:)))./((deltate*2*pi)); 
    figure;colormap(gray);imagesc(flipud(squeeze(abs(one(coil,phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'PHASE  ', 'Coil = ', coil));
    figure;colormap(gray);imagesc(flipud(squeeze((freq(coil,phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'FREQ  ', 'Coil = ', coil));colormap(CMRmap);colorbar
end

one2 = permute(squeeze(mat2gray(abs(one))), [2,3,1]);
one2 = one2*720-360;
freq2 = permute(squeeze(freq), [2,3,1]);
result = cat(4, one2, freq2);
result = rot90(result,1);
%nii= make_nii(result, [2,2,6]);         
%save_nii(nii, ['C:\Users\Ricardo\Desktop\temp\b0map\' measfile(49:55) 'b0map.nii']);
%%

%%PLOT ORIGINAL GRE IMAGES
%%DRAWING LA & PV ROIs FOR OFF-RESONANCE CALCULATIONS
top = zeros(npe, ns);
bottomr= zeros(npe,ns);
for phase = 1:phases %GENERALLY SINGLE PHASE
figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(one(:,phase,:,:))));
    colormap(gray); imagesc(fliplr((finalIm(:,:)))); axis equal; title('TE = 5ms ');

figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,phase,:,:))));
    colormap(gray); imagesc(fliplr((finalIm(:,:)))); axis equal; title('TE = 3.6ms ');hold on 
    
    message1 = sprintf('Left click and hold to begin drawing the LA ROI.');
    uiwait(msgbox(message1));
    [LA, x, y] = roipoly;
    LA = poly2mask(x, y, npe, ns);
    message2 = sprintf('Left click and hold to begin drawing the left PV ROI.');
    uiwait(msgbox(message2));
    [LPV, a, b] = roipoly;
    LPV = poly2mask(a,b,npe,ns);
    message3 = sprintf('Left click and hold to begin drawing the right PV ROI.');
    uiwait(msgbox(message3));
    [RPV, c, d] = roipoly;
    RPV = poly2mask(c,d,npe,ns);

%IF ROIs HAVE ALREADY BEEN DRAWN LOAD HERE AND COMMENT OUT THE ABOVE
% load RPV
% load LA
% load LPV

% %PLOT FREQUENCY AND PHASE MAPS WITH COILS MERGED
% figure;finale(phase,:,:)=mergeCoilImages(squeeze(pha(:,phase,:,:)));
%     colormap(gray); imagesc(fliplr(squeeze(finale(phase,:,:)))); axis equal; title('Phase of Difference ');
% 
% figure;finale(phase,:,:)=mergeCoilImages(squeeze(freq(:, phase,:,:)));
%     colormap(gray); imagesc(fliplr(squeeze(finale(phase,:,:)))); axis equal; title('Freq of Difference ');

    
%USE BERNSTEIN EQUATION TO CALCULATE PHASE
%ENSURE THE NUMBER OF TERMS CORRESPOND TO THE NUMBER OF COILS
    for n = 1:npe
        for k = 1:ns    

if nc == 4            
top = ((abs(one(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(one(2,1,n,k))).^2.*pha(2,1,n,k))+((abs(one(3,1,n,k))).^2.*pha(3,1,n,k))+((abs(one(4,1,n,k))).^2.*pha(4,1,n,k));%+(abs((one(5,1,n,k))).^2.*freq(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
bottom =  (abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2;%;+(abs(one(5,1,n,k))).^2;%+(abs(one(6,1,n,k))).^2;
combo(n,k) = double(top/bottom);
end
if nc == 5
top = ((abs(one(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(one(2,1,n,k))).^2.*pha(2,1,n,k))+((abs(one(3,1,n,k))).^2.*pha(3,1,n,k))+((abs(one(4,1,n,k))).^2.*pha(4,1,n,k))+((abs(one(5,1,n,k))).^2.*pha(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
bottom =  (abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2+(abs(one(5,1,n,k))).^2;%+(abs(one(6,1,n,k))).^2;
combo(n,k) = double(top/bottom);
end
if nc == 6
top = ((abs(one(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(one(2,1,n,k))).^2.*pha(2,1,n,k))+(abs(one(3,1,n,k)).^2.*pha(3,1,n,k))+((abs(one(4,1,n,k))).^2.*pha(4,1,n,k))+(abs((one(5,1,n,k))).^2.*pha(5,1,n,k))+(abs((one(6,1,n,k))).^2.*pha(6,1,n,k));
bottom =  (abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2+(abs(one(5,1,n,k))).^2+(abs(one(6,1,n,k))).^2;
combo(n,k) = double(top/bottom);
end

% up = ((abs(one(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(one(2,1,n,k))).^2.*pha(2,1,n,k))+(abs(one(3,1,n,k)).^2.*pha(3,1,n,k))+((abs(one(4,1,n,k))).^2.*pha(4,1,n,k));%+(abs((one(5,1,n,k))).^2.*pha(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*pha(6,1,n,k));
% down =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(one(5,1,n,k))).^2;%+(abs(one(6,1,n,k))).^2;
% together(n,k) = double(up/down);
% 
% 
% numo = (abs(diff(1,1,n,k))).^2.*angle(diff(1,1,n,k))+(abs(diff(2,1,n,k))).^2.*angle(diff(2,1,n,k))+(abs(diff(3,1,n,k))).^2.*angle(diff(3,1,n,k))+(abs(diff(4,1,n,k))).^2.*angle(diff(4,1,n,k));%+(abs((diff(5,1,n,k))).^2.*;%pha(5,1,n,k));%+(abs((diff(6,1,n,k))).^2.*freq(6,1,n,k));
% deno =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(diff(5,1,n,k))).^2;%+(abs(diff(6,1,n,k))).^2;
% frac(n,k) = double(numo/deno);
        end
    end
phadif =fliplr(combo./(deltate*2*pi)); 


%PLOT RESULTANT FIELD MAP
figure; map = imagesc((phadif)); colormap(gray);  axis equal; title('GRE FIELD MAP - combo '); set(gca,'Visible','off');hold on;
% figure; map = imagesc(fliplr(together)); colormap(gray);  axis equal; title('GRE FIELD MAP - together '); set(gca,'Visible','off');hold on; 
% figure; map = imagesc(fliplr(frac)); colormap(gray);  axis equal; title('GRE FIELD MAP -frac '); set(gca,'Visible','off');hold on;

%OVERLAY THE ROIs ONTO THE PLOTTED FIELD MAP
    plot(x,y, 'r.-', 'MarkerSize', 10)
    plot(a,b, 'b.-', 'MarkerSize', 10)
    plot(c,d, 'g.-', 'MarkerSize', 10)
 
figure; map = imagesc((phadif)); colormap(CMRmap); colorbar; axis equal; title('GRE FIELD MAP - combo '); set(gca,'Visible','off');hold on;  plot(x,y, 'k.-', 'MarkerSize', 10)   
    %GENERATE MASKS IN CORRECT FORMAT
    %CALCULATE MEAN FREQUENCY WITHIN ROI
    LA_filter = double(LA);
%     section = together.*LA_filter;
%     for i = 1:npe; for j = 1:ns; if section(i,j) == 0; section(i,j) = finalIm(i,j); end; end; end
%     imagesc(section)
newlaoffres= sum(sum(LA_filter.*phadif))/sum(sum(LA_filter))
    LA_filter(LA_filter == 0) = NaN;
    leftatrium = phadif.*LA_filter;
    LA_mean = mean(leftatrium(~isnan(leftatrium)));
    
    LPV_filter = double(LPV);
    newlpvoffres= sum(sum(LPV_filter.*phadif))/sum(sum(LPV_filter))
    LPV_filter(LPV_filter == 0) = NaN;
    leftpv = phadif.*LPV_filter;
    LPV_mean = mean(leftpv(~isnan(leftpv)));
    
    
    RPV_filter = double(RPV);
    figure; title('roi checker');
    imagesc((1-RPV_filter).*phadif);
   newrpvoffres= sum(sum(RPV_filter.*phadif))/sum(sum(RPV_filter))
    RPV_filter(RPV_filter == 0) = NaN;
    rightpv = phadif.*RPV_filter;
    
    RPV_mean = mean(rightpv(~isnan(rightpv)));
    
    %TO TAKE VALUES RELATIVE TO LEFT ATRIUM
    offset = 0-LA_mean;
    
    %PRINT OUT THE MEASURED VALUES
    words = "The mean frequency of the LA is: %0.5f (%d) Hz\nThe mean frequency of the left PV is: %0.5f (%0.5f) Hz\nThe mean frequency of the right PV is: %0.5f (%0.5f) Hz";
    sprintf(words,LA_mean,0, LPV_mean, LPV_mean +offset, RPV_mean, RPV_mean + offset)
end


%     LA_filter = double(LA);
% %     section = together.*LA_filter;
% %     for i = 1:npe; for j = 1:ns; if section(i,j) == 0; section(i,j) = finalIm(i,j); end; end; end
% %     imagesc(section)
%     LA_filter(LA_filter == 0) = NaN;
%     leftatrium = squeeze(freq(2,1,:,:)).*LA_filter;
%     LA_mean = mean(leftatrium(~isnan(leftatrium)));
%     
%     LPV_filter = double(LPV);
%     LPV_filter(LPV_filter == 0) = NaN;
%     leftpv = squeeze(freq(2,1,:,:)).*LPV_filter;
%     LPV_mean = mean(leftpv(~isnan(leftpv)));
%     
%     RPV_filter = double(RPV);
%     RPV_filter(RPV_filter == 0) = NaN;
%     rightpv = squeeze(freq(2,1,:,:)).*RPV_filter;
%     RPV_mean = mean(rightpv(~isnan(rightpv)));
%     
%     %TO TAKE VALUES RELATIVE TO LEFT ATRIUM
%     offset = 0-LA_mean;
%     
%     %PRINT OUT THE MEASURED VALUES
%     words = "The mean frequency of the LA is: %0.5f (%d) Hz\nThe mean frequency of the left PV is: %0.5f (%0.5f) Hz\nThe mean frequency of the right PV is: %0.5f (%0.5f) Hz";
%     sprintf(words,LA_mean,0, LPV_mean, LPV_mean +offset, RPV_mean, RPV_mean + offset)

