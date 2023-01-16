
% Dana C. Peters (C) 2015 

% step 1: make corrupted kspace with image phase and delayed kspace lines.
% step 2: estimate delays using Rasche method (to approximately
% anti-parallel projections--measure linear phase.
% step 3: correct kspace delays by adding linear phase 
% Step 4: recon kspace data with iradon to test correction.

%REFERNCES
%Volker Rasche MRM 42 1999 mR fluoroscopy using projection reconstruction
%CB Ahn ZH CHO IEEE  col MI-6 1 1987  a new phase correction method in nmr
%imaging.
%DC Peters, JA Derbyshire, MRM 50 (2003) Centering the projection reconstruction trajectory;
%  

% test by changing roll, delx and dely;  can change Nx and Np.
% known Issues:  

%1. requires removing 0th phase error... ??
%2. Not designed to cope with odd Np yet... easy fix.


% Step 1: Make data with 1) linear phase across image, and 2) delays in kspace in x
% and y.
close all
clear all
Nx=128;
Np=180;
ph=phantom(Nx);
dontusecorr=0; %turn this on to see how uncorrected image looks
roll=200*Nx;% add a phase roll to the image, as might happen with coil phase or shim
for l=1:Nx
    for j=1:Nx
        
        phimag(l,j)=ph(l,j)*exp(-i*2*pi*l/roll) ;
    end
end
imagesc(imag(phimag)); title('imaginary part of test image');  
% look at image with phase roll.
theta = (0:Np/(Np-1): Np)*360/Np;% # of angles for simulation...over 360 deg
[sino,xp] = radon(phimag,theta,Nx);%this is the sinogram data.
delx=3;%delays in kspace...
dely=0;
for p=1:Np
    % the delay depends on the the relative conributions of the gradients.
    % note the abs(cos(theta)...
    del=delx*abs(cos(theta(p)*pi/180))+dely*abs(sin(theta(p)*pi/180));
    for xx=1:Nx
        sino(xx,p)=sino(xx,p)*exp(-i*2*pi*del*xx/Nx);
    end
    totaldel(p)=-del+(Nx/roll)*sin(theta(p)*pi/180); %this will be used to compare to measured later.
end


rawkfft=fftshift(fft(double(sino(:,:)),[],1),1);% fft the radon data into raw kspace sinogram
[Nx,Np]=size(sino);
figure; imagesc(angle(sino));title('original sinogram');
figure, imagesc(abs(rawkfft));title('raw kspace');


%rawkfft is the corrupted data --it is kspace with Nx x Np matrix. the
%number of projections is even and goes through 0 to 2pi.
%STEP 2:   linear phase  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now the processing of the data starts, using the rawkfft as input raw
%kspace sinogram data
linphase=zeros(1,Np);
for p=1:Np
    
    kspalinex(:,p)=(ifft(fftshift(squeeze(rawkfft(:,p))))); %fft each proj into image space
end

for p=1:Np
    
    kspaline=kspalinex(:,p); %fft each proj into image space
    val=0;
    %cho and ahn paper, where the linear phase is estimated
    for m=1:Nx-1
        val=val+kspaline(m)*conj(kspaline(m+1));  %calculate the delta-phase, based on the Ahn and Cho algorithm.
    end
    
    linphase(p)=Nx*angle(val)/(2*pi);  % this is the linear phase is units of pixels of kspace shift
    
end
figure;plot(theta,linphase);% linear phase overall kspace. 
hold; plot(theta, -totaldel, 'g*');
title('linphase'); legend('estimated error', 'theoretical error');hold off;

for p=1:Np/2
    %- gives the object phase
    % + gives the system phase (delay)  must be removed by
    % subtracting...this pixel shift...
    linphasecorr(p)=(linphase(p)+linphase(Np/2+p))/2;% from rasche's paper.
    %prepare for removing the system linear phase from each projection
end
if (dontusecorr)
    linphasecorr=0*linphasecorr;
end
%Step 3: remove delays
jvect=double(transpose([-Nx/2+0.5:1:Nx/2-0.5]/Nx)); % this is right, then 2*pi*pixel shift will work.
figure;plot(theta(1:Np/2), linphasecorr);
title('system linear phase error--in pixel shifts--to be removed'); 
hold off;
clear kspaline
for t=0:Np-1
    kspaline=(fft(fftshift(squeeze(rawkfft(:,t+1))))); % fftd into projection space...
    %if t==2 plot(abs(kspaline));end  for debugging
    
    % projection space now--must not be fftshifted here...
    % fft kspace even and odd projections into projection space.
    num=floor((t/(Np/2))); % get the correct linear phase correction value
    sign=-1; %
    
    phi=(2*pi*[-Nx/2+0.5:1:Nx/2-0.5]/Nx)*(1*sign*linphasecorr(mod(t,Np/2)+1)); %this is the  phase at each point across the image
    
    kspalinexx=kspaline.*(exp(1i*phi')); %perform linear phase correction in projection space
    kspaline2=fftshift(ifft(fftshift(kspalinexx)));  %fft into kspace line...
    phi0=angle(kspaline2(Nx/2+1));%remove 0th order phase--get kspace central phase.
    
    allraw(:,t+1)=exp(-i*phi0)*kspaline2;
    
end

%FINISHED allraw is the correcte kspace!!!

%step 5: test correction by iradon recon..
figure;subplot(1,2,1);
imagesc(real(rawkfft));title('sinogram kspace with delays');
subplot(1,2,2); imagesc(real(allraw));title('corrected raw sinogram space'); %corrected raw kspace...
phi_0=0;

sinocorr=(fftshift(ifft(fftshift(allraw(:,:),1),[],1),1));%fft into sinogram space again.
figure;imagesc(angle(sinocorr));title('corrected sinogram'); 
recon(:,:)=iradon(real(sinocorr),theta, 'linear','Ram-Lak',Np)+1i*iradon(imag(sinocorr),theta, 'linear','Ram-Lak',Np);
% bug in iradon, must perofrm real and imaginary seperately??

figure; imagesc(squeeze(abs(recon(:,:)))); hold on;% corrected image;
title('linear delay corrected image');
axis equal;



