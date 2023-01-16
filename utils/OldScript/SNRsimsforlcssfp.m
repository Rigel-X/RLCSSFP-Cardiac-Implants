% so lets try bssfp first
close all 
clear all


B0=15;
gamma=19.6/3;% read this in De Feyter's Paper.
freq0=gamma*B0; % in MHz Mhz/T
off=0;
freq1=0;%1.3*freq0+off;%lactate
freq2=0;%4.0*freq0+off;%glucose
freq3=0;%5*freq0+off;%water deuterated

T1_1=1800;%blood
T2_1=180;


T1_2=64;%glucose
T2_2=32;
T1_3=320%water
T2_3=12;

TR=4;
TE=TR/2;
TEprep=0;
alpha=40*pi/180;
T1=T1_1;
T2=T2_1;

E1=exp(-TR/T1);
E2=exp(-TR/T2);
phi=0;%offset angle...

j=1;
% different phase offsets
% phaseoffset=0;% this is the 0,0
phaseoffset=180;% this is 0,90
% phaseoffset=2; % this is 0, 180
% phaseoffset=3% this is 0 270
freq1=(1/(TR*.001))*phaseoffset;% shift in frequency, to demonstrate phase offset. 


%for 
for freq=-300:300 %hz
    phi=2*pi*(freq+freq1)*TR*.001;
    i=sqrt(-1);

a=-(1-E1)*E2*sin(alpha);
b=(1-E1)*sin(alpha);
c=E2*(E1-1)*(1+cos(alpha));
d=1-E1*cos(alpha)-(E1-cos(alpha))*E2*E2;

Mssxy(j)=exp(i*(phi/2+pi/2))*(a*exp(-i*phi)+b)./(c*cos(phi)+d);
% Mssxy(j)=(a.*exp(-i*phi)+b)./(c*cos(phi)+d);
Mxyss(j)=((1-E1)*sin(alpha))/(1-(E1-E2)*cos(alpha) -E1*E2);% not used, steady state signal.


freqi(j)=freq;
j=j+1;
end

figure;
plot(freqi,abs(Mssxy),'k-');hold on
%plot(freqi,Mxyss,'g*');

%
M1 = Mssxy;
N=63;
M4(1:N) = M1(252-N:251);
M4(N+1:601) = M1(1:601-N);
N=125;
M2(1:N) = M1(252-N:251);
M2(N+1:601) = M1(1:601-N);
N=188;
M3(1:N) = M1(252-N:251);
M3(N+1:601) = M1(1:601-N);
M1=abs(M1);
M2=abs(M2);
M3=abs(M3);
M4=abs(M4);

figure,plot((M1+M2+M3+M4)/4)
ylim([0,0.2])

for ii = 1:601
    temp(1)=M1(ii);temp(2)=M2(ii);temp(3)=M3(ii);temp(4)=M4(ii);
    mip(ii) = max(temp);
end

figure,plot(mip)
ylim([0,0.2])


%%
close all

figure,plot(freqi,abs(M1),'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.18])

figure,plot(freqi,abs(M4),'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.18])

figure,plot(freqi,abs(M2),'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.18])


figure,plot(freqi,abs(M3),'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.18])


figure,plot(abs((M1+M2+M3+M4))/4,'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.2])

figure,plot(freqi,abs(mip),'white-','LineWidth',2);
grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0.04,0.18])

%%
close all
figure,
for j = 1:4
if j == 1
plot(freqi,abs(M1),'white-','LineWidth',2);
end

if j == 2
plot(freqi,abs(M4),'white-','LineWidth',2);
end

if j == 3
plot(freqi,abs(M2),'white-','LineWidth',2);
end

if j == 4
plot(freqi,abs(M3),'white-','LineWidth',2);
end

grid on
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);
set(gca,'xcol','w','ycol','w')
ylim([0,0.18])

drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'bssfp.gif','gif','loopcount',inf,'Delaytime',0.5)
else
imwrite(I,map,'bssfp.gif','gif','writemode','append','Delaytime',0.5)
end
end