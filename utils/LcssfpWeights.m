% 2022/1/23, jie xiang @yale mrrc


%% Input the Omega map, output four weights map
%[nx,ny]=size(OmegaMap);
% for dyn = 1:4
%     for i = 1:nx
%         for j = 1:ny
%             WeightMap(nx,ny,dyn)=1; % CS-SSFP
%         end
%     end
% end


%% bssfp first
close all 
clear all

% B0=15;
% gamma=19.6/3;% read this in De Feyter's Paper.
% freq0=gamma*B0; % in MHz Mhz/T

T1=320;%water
T2=12;
% T1=1400;%lactate
% T2=200;
% T1=64;%glucose
% T2=32;
TR=10; % ms
freqlim=1000/TR; %Hz
alpha=18*pi/180; % optimal flip angle = arccos(exp(-TR/T1))

% different phase offsets
phaseoffset = [0,0.25,0.5,0.75]; % four dynamics
% phaseoffset = [0,1/3,2/3]; % three dynamics
% phaseoffset = [0,1/2]; % two dynamics
% phaseoffset = [0]; % one dynamic
measure = size(phaseoffset,2);
StepHz = 0.01;

E1=exp(-TR/T1);
E2=exp(-TR/T2);
a=-(1-E1)*E2*sin(alpha);
b=(1-E1)*sin(alpha);
c=E2*(E1-1)*(1+cos(alpha));
d=1-E1*cos(alpha)-(E1-cos(alpha))*E2*E2;

freq1=(1/(TR*.001)).*phaseoffset;% shift in frequency, to demonstrate phase offset. 
j=1;

for freq=-freqlim:StepHz:freqlim %hz
    phi=2*pi*(freq+freq1)*TR*.001;
    Mssxy(j,:)=(a.*exp(-i*phi)+b)./(c*cos(phi)+d);
    freqi(j)=freq;
    j=j+1;
end

Wgts=zeros(2*freqlim/StepHz+1,measure);
[~,minpos] = min(Mssxy(freqlim/StepHz+2:2*freqlim/StepHz+1,measure));
%MinPos = round(minpos*measure/4);
MinPos = round(minpos*measure/2-minpos);
TotLen = measure*minpos;
WtsLen = 2*minpos;
Nloop = 0;
flag=0;
while flag == 0
    for jj = 1:minpos
        % left increasing
        if freqlim/StepHz+1+MinPos+jj+measure*Nloop*TotLen > 2*freqlim/StepHz+1
            flag=1;
            break
        end
        Wgts(round(freqlim/StepHz+MinPos+jj+measure*Nloop*TotLen)+1,1)=jj/minpos;
        Wgts(round(freqlim/StepHz-MinPos-jj-measure*Nloop*TotLen)+1,1)=jj/minpos;
        % right decreasing
        if freqlim/StepHz+1+MinPos+minpos+jj+measure*Nloop*TotLen > 2*freqlim/StepHz+1
            flag=1;
            break
        end        
        Wgts(round(freqlim/StepHz+MinPos+minpos+jj+measure*Nloop*TotLen)+1,1)=1-jj/minpos;
        Wgts(round(freqlim/StepHz-MinPos-minpos-jj-measure*Nloop*TotLen)+1,1)=1-jj/minpos;
    end
    Nloop = Nloop+1;
end

for j=2:measure
    Wgts(1:2*freqlim/StepHz+1-(j-1)*minpos,j) = Wgts((j-1)*minpos+1:2*freqlim/StepHz+1,1);
    Wgts(2*freqlim/StepHz+2-(j-1)*minpos:2*freqlim/StepHz+1,j) = Wgts(2:(j-1)*minpos+1,1);
end

% Wgts(1:2*freqlim+1-MinPos,2) = Wgts(MinPos+1:2*freqlim+1,1);
% Wgts(2*freqlim+2-MinPos:2*freqlim+1,2) = Wgts(2:MinPos+1,1);
% 
% Wgts(1:2*freqlim+1-2*MinPos,3) = Wgts(2*MinPos+1:2*freqlim+1,1);
% Wgts(2*freqlim+2-2*MinPos:2*freqlim+1,3) = Wgts(2:2*MinPos+1,1);
% 
% Wgts(1:2*freqlim+1-3*MinPos,4) = Wgts(3*MinPos+1:2*freqlim+1,1);
% Wgts(2*freqlim+2-3*MinPos:2*freqlim+1,4) = Wgts(2:3*MinPos+1,1);

figure,subplot(2,1,1)
plot(freqi,abs(Mssxy(:,1)),'r-');hold on
plot(freqi,Wgts(:,1)/50+0.11,'r--')
if measure>1
plot(freqi,abs(Mssxy(:,2)),'b-');
plot(freqi,Wgts(:,2)/50+0.11,'b--')
end
if measure>2
plot(freqi,abs(Mssxy(:,3)),'y-');
plot(freqi,Wgts(:,3)/50+0.11,'y--')
end
if measure>3
plot(freqi,abs(Mssxy(:,4)),'g-');%legend('0','90','180','270')
plot(freqi,Wgts(:,4)/50+0.11,'g--')
end

%% compare linear combine and maximum intensity projection
subplot(2,1,2)
WtedM=Wgts.*Mssxy;

if measure == 4
M1=squeeze(abs(Mssxy(:,1)));M2=squeeze(abs(Mssxy(:,2)));M3=squeeze(abs(Mssxy(:,3)));M4=squeeze(abs(Mssxy(:,4)));
linearc = (M1+M2+M3+M4)/4; 
sosc = sqrt(M1.^2+M2.^2+M3.^2+M4.^2)/2.2; 
for ii = 1:2*freqlim/StepHz+1
    temp(1)=M1(ii);temp(2)=M2(ii);temp(3)=M3(ii);temp(4)=M4(ii);
    mip(ii) = max(temp);
end
M1=squeeze(abs(WtedM(:,1)));M2=squeeze(abs(WtedM(:,2)));M3=squeeze(abs(WtedM(:,3)));M4=squeeze(abs(WtedM(:,4)));
fieldc = M1+M2+M3+M4;
end

if measure == 3
M1=squeeze(abs(Mssxy(:,1)));M2=squeeze(abs(Mssxy(:,2)));M3=squeeze(abs(Mssxy(:,3)));
linearc = (M1+M2+M3)/2.5; 
sosc = sqrt(M1.^2+M2.^2+M3.^2)/2.2; 
for ii = 1:2*freqlim/StepHz+1
temp(1)=M1(ii);temp(2)=M2(ii);temp(3)=M3(ii);
mip(ii) = max(temp);
end
M1=squeeze(abs(WtedM(:,1)));M2=squeeze(abs(WtedM(:,2)));M3=squeeze(abs(WtedM(:,3)));
fieldc = M1+M2+M3;
end

if measure == 2
M1=squeeze(abs(Mssxy(:,1)));M2=squeeze(abs(Mssxy(:,2)));
linearc = (M1+M2)/4; 
sosc = sqrt(M1.^2+M2.^2)/2; 
for ii = 1:2*freqlim/StepHz+1
    temp(1)=M1(ii);temp(2)=M2(ii);
    mip(ii) = max(temp);
end
M1=squeeze(abs(WtedM(:,1)));M2=squeeze(abs(WtedM(:,2)));
fieldc = M1+M2;
end

if measure == 1
M1=squeeze(abs(Mssxy(:,1)));
linearc = M1; 
sosc = sqrt(M1.^2); 
for ii = 1:2*freqlim/StepHz+1
    temp(1)=M1(ii);
    mip(ii) = max(temp);
end
M1=squeeze(abs(WtedM(:,1)));
fieldc = M1;
end


DelOverMaxl = (max(linearc)-min(linearc))/max(linearc);
DelOverMaxs = (max(sosc)-min(sosc))/max(sosc);
DelOverMaxm = (max(mip)-min(mip))/max(mip);
DelOverMaxf = (max(fieldc)-min(fieldc))/max(fieldc);
plot(freqi,linearc,'b'),hold on
plot(freqi,sosc,'m'),
plot(freqi,mip/1.2,'g'),hold on
plot(freqi,fieldc,'r'),ylim([0,0.2]);

% figure,plot(freqi,(M2+M3)/2),ylim([0,0.2]), title('linear combine')
legend(['linear combine, Del over Max = ', num2str(DelOverMaxl)],['sos combine, Del over Max = ', num2str(DelOverMaxs)],['Max intensity proj, Del over Max = ', num2str(DelOverMaxm)],['field map weighted combine, Del over Max = ', num2str(DelOverMaxf)])
xlabel('frequency, Hz')


