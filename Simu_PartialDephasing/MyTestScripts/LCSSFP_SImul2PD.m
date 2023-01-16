%% bssfp response profile
clear
close all
load('JXPD.mat')
T1=1000;%water
T2=150;
% T1=1400;%lactate
% T2=200;
% T1=64;%glucose
% T2=32;
TR=3.5; % ms
freqlim=2100/TR; %Hz
alpha=35*pi/180; % optimal flip angle = arccos(exp(-TR/T1))

% different phase offsets
phaseoffset = [0,0.25,0.5,0.75]; % four dynamics
% phaseoffset = [0,1/3,2/3]; % three dynamics
% phaseoffset = [0,1/2]; % two dynamics
% phaseoffset = [0]; % one dynamic
StepHz = 0.01;

E1=exp(-TR/T1);
E2=exp(-TR/T2);
a=-(1-E1)*E2*sin(alpha);
b=(1-E1)*sin(alpha);
c=E2*(E1-1)*(1+cos(alpha));
d=1-E1*cos(alpha)-(E1-cos(alpha))*E2*E2;
freq1=(1/(TR*.001)).*phaseoffset;% shift in frequency, to demonstrate phase offset. 
j=1;
phasecycleNumber = 1;
freqi = zeros(1,1+2*freqlim);
for freq=-freqlim:StepHz:freqlim %hz
    phi=2*pi*(freq+freq1(phasecycleNumber))*TR*.001;
    Mssxy(j,:)=(a.*exp(-1i*phi)+b)./(c*cos(phi)+d);
    freqi(j)=freq;
    j=j+1;
end
min(abs(Mssxy(:)))
max(abs(Mssxy(:)))
figure,plot(freqi,abs(Mssxy),'r-');
figure,plot(freqi,angle(Mssxy),'r-');

% generate intensity mask
BandMask0 = zeros(4,192,192);
BandMask1 = zeros(4,192,192);
BandMask2 = zeros(4,192,192);
Fmap = zeros(192,192);

npc = size(phaseoffset,2);
for i = 1:192
    for j = 1:192
        for Npc = 1:npc           
            dist = sqrt(i^2+j^2);
            Fmap(i,j) = 1000*dist^(1/7);
            phi=2*pi*(Fmap(i,j)+freq1(Npc))*TR*.001;
            BandMask0(Npc,i,j)=(a.*exp(-1i*phi)+b)./(c*cos(phi)+d);
            
            phieff = mod(phi,4*pi)/(4*pi)*720;
            [~,Index] = min(abs(phieff-xq));%Index
            BandMask1(Npc,i,j)=PD00(Index);
            BandMask2(Npc,i,j)=PD600(Index);
            
            BandMask(Npc,i,j) = BandMask1(Npc,i,j);
            
        end
        test(i,j) = max(max(abs(BandMask(:,i,j))));
    end
end

figure,imshow(Fmap,[])
min(abs(BandMask(:)))
max(abs(BandMask(:)))
figure,for i = 1:npc; subplot(2,2,i),imshow(squeeze(BandMask(i,:,:)),[0,max(abs(BandMask(:)))]);end
test2 = 0;
for i = 1:npc; test2 = test2 + BandMask(i,:,:); end
figure,subplot(1,2,1),imshow(squeeze(test),[0,4*max(abs(BandMask(:)))]),title('MIP');subplot(1,2,2),imshow(squeeze(abs(test2)),[0,4*max(abs(BandMask(:)))]),title('lc');

phtmask = (imresize(phantom,192/256)~=0);

figure,
for i = 1:npc
    PhtPacemaker(i,:,:) = imresize(phantom,192/256).*squeeze(BandMask(i,:,:)).*phtmask;
    subplot(2,2,i),imshow(abs(squeeze(PhtPacemaker(i,:,:))),[]);
end

if npc == 4
figure,imshow(abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:)))+abs(squeeze(PhtPacemaker(3,:,:)))+abs(squeeze(PhtPacemaker(4,:,:))),[])
elseif npc == 3
figure,imshow(abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:)))+abs(squeeze(PhtPacemaker(3,:,:))),[])
elseif npc == 2
figure,imshow(abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:))),[])
end


%%
% PhtPacemaker = PhtPacemakerPD60;
if npc == 4
finalImg = abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:)))+abs(squeeze(PhtPacemaker(3,:,:)))+abs(squeeze(PhtPacemaker(4,:,:)));
elseif  npc == 3
finalImg = abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:)))+abs(squeeze(PhtPacemaker(3,:,:)));
elseif  npc == 2
finalImg = abs(squeeze(PhtPacemaker(1,:,:)))+abs(squeeze(PhtPacemaker(2,:,:)));
end

% finalImg = abs(squeeze(PhtPacemaker(1,:,:))+(squeeze(PhtPacemaker(2,:,:)))+(squeeze(PhtPacemaker(3,:,:)))+(squeeze(PhtPacemaker(4,:,:))));
% finalImg = abs(squeeze(PhtPacemaker(1,:,:)));
load('F:\CVCoding\gridding_lcssfp\Simulation\ksptest48.mat', 'Refimg')
Ref = rot90(squeeze(abs(Refimg(:,:,10))),2);
Ref(imresize(phantom,192/256)==0)=1;
Tru = imresize(phantom,192/256);
Tru(imresize(phantom,192/256)==0)=1;
finalImg = finalImg.*abs(Ref)./abs(Tru);
figure,imshow(finalImg,[])

%% generate phantom and ksp
x = 2.5:2.5:360;
xq = 0.01:0.01:360;
PD00 = interp1(x,PD0V0,xq);
PD00(1:249)=PD00(250);

temp(18001:54000)=PD00;
figure,plot(xq,abs(PD00))
for k = 1:18000
temp(k) = abs(PD00(18001-k)).*exp(1i*(pi-angle(PD00(18001-k))));
end
for k = 54001:72000
temp(k) = abs(PD00(k-54000)).*exp(1i*(-pi-angle(PD00(90001-k))));
end
PD00 = temp;
figure,plot(abs(temp))
figure,plot(angle(temp))

PD600 = interp1(x,PD60V0,xq);
PD600(1:249)=PD600(250);
temp(18001:54000)=PD600;
figure,plot(xq,abs(PD600))
for k = 1:18000
temp(k) = abs(PD600(18001-k)).*exp(1i*(pi-angle(PD600(18001-k))));
end
for k = 54001:72000
temp(k) = abs(PD600(k-54000)).*exp(1i*(-pi-angle(PD600(90001-k))));
end
PD600 = temp;
xq = 0.01:0.01:720;
save('JXPD.mat','PD0V0','PD60V0','PD00','PD600','PhtPacemakerPD0','PhtPacemakerPD60','x','xq')


%%
temp = angle(Mssxy(60001:117100,1))';
xx = 0.00630472854:0.00630472854:360;
temp = interp1(xx,temp,xq);
figure,plot(temp)
for k = 1:72000
dddd(k) = abs(PD00(k)).*exp(1i*temp(k));
end
dddd(1:99)=PD00(100);
PD00 = dddd;
for k = 1:72000
PD600(k) = abs(PD600(k)).*exp(1i*temp(k));
end
PD600(1:99)=PD600(100);
figure,plot(abs(PD00))
figure,plot(angle(PD00))
save('JXPD.mat','PD0V0','PD60V0','PD00','PD600','PhtPacemakerPD0','PhtPacemakerPD60','x','xq')


