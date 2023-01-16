clear all

%LOAD IN IMG/NIFTI FILE
nii=load_nii('C:\Users\Ricardo\Desktop\JaimesFiles\Phantoms\0812_bottlephantoms\48P\LC_SSFP_STD_NOINT_48_FOUR.img');
image = nii.img;

%DRAW ROI
% artefact = roipoly;
% signal = roipoly;

%LOAD ROI IF ALREADY SELECTED 
load artefact
load signalphan

%%CONVERT MASK TO SUITABLE FORMAT
sig_fil = double(signal);
art_fil = double(artefact);
sig_fil(sig_fil == 0) = NaN;
art_fil(art_fil == 0) = NaN;

%SCALE IMAGE FROM 0 TO 1
scaale = image/(max(max(max(image))));

% show the image
for idx = 1:7
    figure; imagesc(image(:,:,idx));colormap(gray);
end

%SELECT ROI FROM SCALED IMAGE
bottle = scaale.*sig_fil;
streak = scaale.*art_fil;

%CALCULATE MEASUREMENTS
bottle_mean1 = mean(bottle(~isnan(bottle)));
streak_mean1 = mean(streak(~isnan(streak)));
bottle_max = max(bottle(~isnan(bottle)));
streak_max = max(streak(~isnan(streak)));
bottle_min = min(bottle(~isnan(bottle)));
streak_min = min(streak(~isnan(streak)));
bottle_total = sum(bottle(~isnan(bottle)));
streak_total = sum(streak(~isnan(streak)));
bottle_sd = std(bottle(~isnan(bottle)));
streak_sd = std(streak(~isnan(streak)));


%PRINT OUT MEASURED VALUES
words = "The max signal is: %0.5f  \nThe min signal is: %0.5f \nThe mean signal  is: %0.5f \nThe standard deviation of signal is: %0.5f \nThe total signal is: %0.5f ";
words2 ="The max artefact is: %0.5f \nThe min artefact is: %0.5f \nThe mean artefact  is: %0.5f  \nThe standard deviation of artefact is: %0.5f \nThe total artefact is: %0.5f";
sprintf(words,abs(bottle_max),abs(bottle_min),abs(bottle_mean1),abs(bottle_sd),abs(bottle_total))
sprintf(words2, abs(streak_max), abs(streak_min),abs(streak_mean1), abs(streak_sd), abs(streak_total))
