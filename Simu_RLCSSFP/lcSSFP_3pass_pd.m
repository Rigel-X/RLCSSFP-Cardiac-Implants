% 2022/10/9, simulate 3-pass lc-ssfp with partial dephasing(PD) to show
% that PD help reduce flow artifacts as well as banding artifacts
% superimposed bssfp banding mask with PD
% jie xiang @yale mrrc

clear all
close all
clc

load('lcssfp_simudata.mat','CineImg')

%% bssfp response profile
T1=1800;T2=200;%water
TR=3.5; % ms
freqlim=2100/TR; %Hz
alpha=35*pi/180; % optimal flip angle = arccos(exp(-TR/T1))

% different phase offsets
% phaseoffset = [0,0.25,0.5,0.75]; % four dynamics
phaseoffset = [0,1/3,2/3]; % three dynamics
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

%% create bSSFP intensity mask
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
CinePmk = zeros(Nx,Ny,NoPhase,3);
for i = 1:Nx
    for j = 1:Ny
        for Npc = 1:3           
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
figure,for i = 1:3; subplot(2,2,i),imshow(squeeze(BandMask(i,:,:)),[0,max(abs(BandMask(:)))]);end
test2 = BandMask(1,:,:)+BandMask(2,:,:)+BandMask(3,:,:);
figure,subplot(1,2,1),imshow(squeeze(test),[0,4*max(abs(BandMask(:)))]),title('MIP');subplot(1,2,2),imshow(squeeze(abs(test2)),[0,4*max(abs(BandMask(:)))]),title('lc');

figure,set(gcf,'Color',[1 1 1]);
for j = 1:30
    subplot(2,2,1),imshow(brighten(cat(2,abs(squeeze(CinePmk(:,:,j,1))),abs(squeeze(CinePmk(:,:,j,2))),abs(squeeze(CinePmk(:,:,j,3)))),0.4),[]),title(['one phase cycle, phase = ',num2str(j)])
    subplot(2,2,4),imshow(brighten(1*abs(squeeze(CinePmk(:,:,j,1))+squeeze(CinePmk(:,:,j,2))+squeeze(CinePmk(:,:,j,3)))+abs(squeeze(CinePmk(:,:,j,1)))+abs(squeeze(CinePmk(:,:,j,2)))+abs(squeeze(CinePmk(:,:,j,3))),0.4),[]),title('lcSSFP'),pause(0.05)
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
save('lcssfp_pd_simudata.mat','CinePmk','CineImg','-append')

%% PD response
PDFac = ((max(abs(Mssxy))-min(abs(Mssxy)))/(max(abs(PD600))-min(abs(PD600))));
Mssxy_pd60 = ((abs(PD600)-min(abs(PD600)))*PDFac + min(abs(Mssxy))).*exp(1i*angle(PD600));
figure,subplot(2,2,1),plot(abs(Mssxy_pd60)),subplot(2,2,2),plot(angle(Mssxy_pd60))
save('PD_bSSFP_res.mat','Mssxy_pd60','-append')


%% create bSSFP_PD intensity mask
% banding mask
Nx = 192; Ny = 192;NoPhase = 30;
BandMask_pd = zeros(4,Nx,Ny);
Fmap = zeros(Nx,Ny);
CinePmk_pd = zeros(Nx,Ny,NoPhase,3);
for i = 1:Nx
    for j = 1:Ny
        for Npc = 1:3           
            dist = ((i-PmkX)^2+(j-PmkY)^2);
            Fmap(i,j) = 600*dist^(1/7);
            FreLoc = max(floor(mod(Fmap(i,j)+(Npc-1)*(571.4*2/3),571)/571.4*72000),1);
            BandMask_pd(Npc,i,j)=Mssxy_pd60(FreLoc);
            for nphases = 1:NoPhase
                CinePmk_pd(i,j,nphases,Npc) = abs(CineImg(i,j,nphases)).*abs(BandMask_pd(Npc,i,j)).*exp(1i*angle(BandMask_pd(Npc,i,j))).*DephMask(i,j);
            end
        end
        test(i,j) = max(max(abs(BandMask_pd(:,i,j))));
    end
end
figure,imshow(Fmap,[])
min(abs(BandMask_pd(:)))
max(abs(BandMask_pd(:)))
figure,for i = 1:3; subplot(2,2,i),imshow(squeeze(BandMask_pd(i,:,:)),[0,max(abs(BandMask_pd(:)))]);end
test2 = BandMask_pd(1,:,:)+BandMask_pd(2,:,:)+BandMask_pd(3,:,:);
figure,subplot(1,2,1),imshow(squeeze(test),[0,4*max(abs(BandMask_pd(:)))]),title('MIP');subplot(1,2,2),imshow(squeeze(abs(test2)),[0,4*max(abs(BandMask_pd(:)))]),title('lc');

figure,set(gcf,'Color',[1 1 1]);
for j = 1:30
    subplot(2,2,1),imshow(brighten(cat(2,abs(squeeze(CinePmk_pd(:,:,j,1))),abs(squeeze(CinePmk_pd(:,:,j,2))),abs(squeeze(CinePmk_pd(:,:,j,3)))),0.4),[]),title(['one phase cycle, phase = ',num2str(j)])
    subplot(2,2,4),imshow(brighten(0.5*abs(squeeze(CinePmk(:,:,j,1))+squeeze(CinePmk(:,:,j,2))+squeeze(CinePmk(:,:,j,3)))+abs(squeeze(CinePmk_pd(:,:,j,1)))+abs(squeeze(CinePmk_pd(:,:,j,2)))+abs(squeeze(CinePmk_pd(:,:,j,3))),0.4),[]),title('lcSSFP'),pause(0.05)
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
save('lcssfp_pd_simudata.mat','CinePmk_pd','-append')



