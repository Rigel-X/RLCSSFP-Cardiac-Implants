% gridding based on hargreaves spiral gridding code.
clear all;

correct_linear=1;
correct_zero=0;
remove_allphase=0;
if (remove_allphase==1)
    correct_linear=1;
end

%set initial parameters
res=160;    %image resolution
kwidth=3;   %kernel width of convolution windo
oversmpl=1.5;    %oversampling-factor gridsize/FOV/resolution
path(path,'C:\Users\Ricardo\Desktop\0914_YR2twix');
folder = 'C:\Users\Ricardo\Desktop\0914_YR2twix';
filepath = '';
measfile = 'C:\Users\Ricardo\Desktop\0914_YR2twix\meas_MID118_LC_Rad_180_FID28362.dat'
dispopt = 'off';
filename = 'meas_MID118_LC_Rad_180_FID28362.dat';   % radial 128 spoke, long TR, 100 mm FOV
%change file name
[rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17phasenew(measfile, dispopt);
%rawdata= ReadSiemensMeasVB17([filepath folder filename],'off');

phases = sMDH.sLC.ushPhase;
nc=sMDH.sLC.ulChannelId; %sMDH.ushUsedChannels;
for cc=1:nc, 
    data(:,:,cc) = (rawdata(cc:nc:end-nc,:));  
end
nc=sMDH.sLC.ulChannelId; %sMDH.ushUsedChannels;
ksamps=permute(data,[3 2 1]);
%ksamps=ksamps(:,:,1:end/2);
coils = size(ksamps,1)
Nx = size(ksamps,2)
Np = size(ksamps,3)/phases;
ksampsnew=ksamps(:,:,1:Np);
phival=[0:pi/Np:pi-pi/Np];
% phival=atan2(gradient_strength_mT_m(:,2),  gradient_strength_mT_m(:,1)); 
nunder=1;% undersampling factor
size(ksampsnew);
ksampsnew= ksampsnew(:,:,1:nunder:end);
phival= phival(1:nunder:end);


phases=1;% using single phase data now.


% density compensation filter...
%create Shepp-Logan dcf (looks smoother than ramlak)
k=[-.5:1/Nx:.5-1/Nx];
dcf1=abs(k).*sinc(k);
dcf=repmat(dcf1,[Np 1]); 

%create dcf-matrix (Ram-lak-filter) (same for all coildata)
% dcf1=[.5:-.5/(Nx/2):0,0.5/(Nx/2):0.5/(Nx/2):0.5-0.5/(Nx/2)];
% % dcf1=abs(fftshift(fft(ramLak_mex(Nx))));
% dcf2=dcf1/max(dcf1);
% dcf=repmat(dcf2,[Np 1]); %384*39 matrix
% figure;plot(dcf(1,:));

%end % if correct_all
%locate kspace samples in complex space (same for all coildata)
loc=zeros(Nx,Np);

for j=[1:Nx] 
    for k=[1:Np] 
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
       
        pha=phival(k);
        loc(j,k)=mag*exp(i*pha);
      %  plot(real(loc),imag(loc),'g*');
    end
end

figure;
%perform gridding (different for all coildata) and create image with
%fourier transform
dat=zeros(phases,coils,res,res); %all k-space data (empty matrix)
im=zeros(phases,coils,res,res);  %all image data (empty matrix)
finalIm=zeros(phases,res,res);   %combined image data from coils (empty matrix)
for ph=1:phases
    for cl=1:coils
        
        size(squeeze(ksampsnew(cl,:,:)))
        dat(ph,cl,:,:) = gridkb(loc,squeeze(ksampsnew(cl,:,:)),dcf',res,kwidth,oversmpl);
        im(ph,cl,:,:) = (fft2(squeeze(dat(ph,cl,:,:))));
        %now "dat" contains all gridded k-space data, "im" contains all
        %image data
    end
    finalIm(ph,:,:)=mergeCoilImages(squeeze(abs(im(ph,:,:,:))));
    colormap(gray); imagesc(fftshift(squeeze(finalIm(ph,:,:))));
   
end 
