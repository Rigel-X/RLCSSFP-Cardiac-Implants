% 2022/1/22, phase correction of radial lc-ssfp
% jie xiang @yale mrrc

function [KspaCored] = LcssfpPhaseCor(KspaToCor)
%remove the linear and 0-order phase
%KspaToCor=allraw;
[Nx,Np]=size(KspaToCor);
linphase=zeros(1,Np);
kspalinep=zeros(Np,Nx);
valflag=zeros(Np,Nx-1);
%% Calculate delay
for p=1:Np
    kspalinep(p,:)=(ifft(fftshift(squeeze(KspaToCor(:,p))))); %fft each proj into image space
    val=0;
    for m=1:Nx-1
        valflag(p,m) = kspalinep(p,m)*conj(kspalinep(p,m+1));
        val=val+valflag(p,m);  %calculate the delta-phase, based on the Ahn and Cho algorithm.
    end    
    linphase(p)=Nx*angle(val)/(2*pi);  % this is the linear phase in units of pixels of kspace shift   
end

linphasecorr=zeros(1,Np);
for p=1:(Np-1)/2
    %- gives the object phase, off-resonance related linear phase error
    % + gives the system phase (delay), system imperfection, not rely on
    % its polarity, must be removed by subtracting...this pixel shift...
    linphasecorr(p)=(linphase(p)+linphase((Np+1)/2+p))/2;% from rasche's paper.
    %prepare for removing the system linear phase from each projection
end

plot(linphasecorr)

%% remove delays
KspaCored = zeros(Nx,Np);
for t=0:Np-1
    %DCP added fftshift here....
    kspalined=fftshift(fft(fftshift(squeeze(KspaToCor(:,t+1))))); % fftd into projection space...
    % kspalined=fft(fftshift(squeeze(KspaToCor(:,t+1)))); % fftd into projection space...
    % projection space now--must not be fftshifted here...
    % fft kspace even and odd projections into projection space.
    % num=floor((t/(Np/2))); % get the correct linear phase correction value
    sign=-1; %
    if t== (Np-1)/2
        a=1+1;
    else
        phi=(2*pi*[-Nx/2+0.5:1:Nx/2-0.5]/Nx)*(1*sign*linphasecorr(mod(t,(Np+1)/2)+1)); %this is the  phase at each point across the image
    
    kspalinexx=kspalined.*(exp(1i*phi')); %perform linear phase correction in projection space
    kspaline2=fftshift(ifft(fftshift(kspalinexx)));  %fft into kspace
 
  
%  kspaline2=(ifft((kspalinexx)));  %fft into kspace line...
    phi0=angle(kspaline2(Nx/2+1)); %remove 0th order phase--get kspace central phase.
    end
    %phi0=0;
    KspaCored(:,t+1)=exp(-1i*phi0)*kspaline2;
    
    
end
%%
% figure,
% subplot(1,2,1),imshow(abs(KspaToCor),[]),title('kspace to be corrected');
% subplot(1,2,2),imshow(abs(KspaCored),[]),title('corrected kspace')
% close all
% figure,plot(abs(KspaToCor(:,Np))),title('kspatocor');
% figure,plot(linphase),title('linphase');
% figure,plot(linphasecorr),title('linphasecorr');
% figure,plot(abs(KspaCored(:,Np))),title('kspacored');
% 
% for i =1:Np
% [~,MaxP(i)]=max(abs(KspaCored(:,i)));
% end
% for i =1:Np
% [~,MaxP2(i)]=max(abs(KspaToCor(:,i)));
% end
% figure,plot(MaxP),hold on
% plot(MaxP2),title('maximum position')

end
