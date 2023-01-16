
folder = 'C:\Users\Dana Peters\Documents\';
filepath = '';

%GRE MAP FILES

% measfile = 'C:\Users\Ricardo\Desktop\1005_YR3\meas_MID49_gre3_6te_FID31657.dat'

dispopt = 'on';
[rawdata,newraw, loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17phasenew(measfile, dispopt);
phases = sMDH.sLC.ushPhase;         %phases
nc=sMDH.sLC.ulChannelId;            %coils
npe = sMDH.sLC.ushLine;             %lines
ns = sMDH.ushSamplesInScan;         %samples
deltate = (5-3.6)/1000; %seconds

% one = zeros(nc, phases, npe, ns);
two = zeros(nc, phases, npe, ns);
diff = zeros(nc, phases, npe, ns);
freq = zeros(nc, phases, npe, ns);
pha = zeros(nc, phases, npe, ns);

for coil = 1:nc
    for phase = 1:phases
%         one(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(newraw(coil, phase, :, :)))));
        two(coil, phase,:,:) = ifftshift(ifft2(ifftshift(squeeze(newraw(coil, phase, :, :)))));
       
     
    end
end
k = 1
% phase =1;
for coil = 1:nc
    diff(coil,1,:,:) = (one(coil,1,:,:)) - (two(coil,1,:,:));
%     phasediff = angle(one(coil,phase,:,:)) - angle(two(coil,phase,:,:));
%     figure;colormap(gray);imagesc(squeeze(angle(diff(coil,phase,:,:)))); axis equal; title(sprintf( "%s%s%d", 'DIFF  ', 'Coil = ', coil));
%     figure;imagesc(squeeze(phasediff(coil,phase,:,:))); axis equal; title(sprintf( "%s%s%d", 'PhaseDIFF  ', 'Coil = ', nc));
    
    pha(coil,1,:,:) = angle((diff(coil,1,:,:)));%sqrt(-1)*
end
% for coil = 1:nc
% while k < 3
% for i = 2:npe
% for j = 2:ns
% 
% chnge = pha(coil,phase,i,j)-pha(coil,phase,i,j-1);
%     if abs(chnge) > pi && chnge>0
%         pha(coil,phase,i,j) = pha(coil, phase, i, j) - (pi);
%         newchnge = pha(coil,phase,i,j) - pha(coil,phase, i,j-1);
% %     if abs(newchnge) > 1 
% %         pha(coil,phase,i,j) = pha(coil,phase,i,j) -(pi);
% %     end
%     end
%     if abs(chnge) > pi && chnge<0
%         pha(coil,phase,i,j) = pha(coil,phase,i,j) +(pi);
%         newchnge = pha(coil,phase,i,j) - pha(coil,phase, i,j-1);
% %     if abs(newchnge) > 1 
% %         pha(coil,phase,i,j) = pha(coil,phase,i,j) +(pi);
% %     end
%     end
% 
% 
% 
% chnge = pha(coil,phase,i,j)-pha(coil,phase,i-1,j);
%     if abs(chnge) > pi && chnge>0
%         pha(coil,phase,i,j) = pha(coil, phase, i, j) - (pi);
%         newchnge = pha(coil,phase,i,j) - pha(coil,phase, i-1,j);
% %     if newchnge > pi 
% %         pha(coil,phase,i,j) = pha(coil,phase,i,j) -(pi);
% %     end
%     end
%     if abs(chnge) > pi && chnge<0
%         pha(coil,phase,i,j) = pha(coil,phase,i,j) + (pi);
%         newchnge = pha(coil,phase,i,j) - pha(coil,phase, i-1,j);
% %     if newchnge > pi 
% %         pha(coil,phase,i,j) = pha(coil,phase,i,j) +(pi);
% %     end
% 
%     end
% end 
% end
% k=k+1
% end
% end


% wrap the image
% wimg = wrapTo2Pi(pha(1,1,:,:));

% tic;
% unwrap_img = unwrap_phase(pha(1,1,:,:));
% toc;
for coil = 1:nc
    freq(coil,phase,:,:) = ((pha(coil,phase,:,:)))./((deltate*2*pi));
    
    figure;colormap(gray);imagesc(flipud(squeeze(abs(pha(coil,phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'PHASE  ', 'Coil = ', coil));
%     figure;colormap(gray);imagesc(fliplr(squeeze(phac(1,1,:,:)))); axis equal;
    figure;colormap(gray);imagesc(flipud(squeeze((freq(coil,phase,:,:))))); axis equal; title(sprintf( "%s%s%d", 'FREQ  ', 'Coil = ', coil));
end


diff=permute(diff,[2 1 3 4]);
% freq=permute(freq,[2 1 3 4]);
top = zeros(npe, ns);
bottomr= zeros(npe,ns);
for phase = 1:phases
figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(one(:,phase,:,:))));
    colormap(gray); imagesc(flipud((finalIm(:,:)))); axis equal; title('TE = 5ms ');

figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,phase,:,:))));
    colormap(gray); imagesc(flipud((finalIm(:,:)))); axis equal; title('TE = 3.6ms ');hold on 
%     message1 = sprintf('Left click and hold to begin drawing the LA ROI.');
%     uiwait(msgbox(message1));
%     [LA, x, y] = roipoly;
%     LA = poly2mask(x, y, npe, ns);
%     message2 = sprintf('Left click and hold to begin drawing the left PV ROI.');
%     uiwait(msgbox(message2));
%     [LPV, a, b] = roipoly;
%     LPV = poly2mask(a,b,npe,ns);
%     message3 = sprintf('Left click and hold to begin drawing the right PV ROI.');
%     uiwait(msgbox(message3));
%     [RPV, c, d] = roipoly;
%     RPV = poly2mask(c,d,npe,ns);
% %     
    for n = 1:npe
        for k = 1:ns
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(angle(diff(phase,:,:,:))));
%     colormap(gray); imagesc(fliplr(squeeze(finalIm(phase,:,:)))); axis equal; title('Phase of Difference ');
% 
% figure;finalIm(phase,:,:)=mergeCoilImages(squeeze(freq(phase,:,:,:)));
%     colormap(gray); imagesc(fliplr(squeeze(finalIm(phase,:,:)))); axis equal; title('Freq of Difference ');
% combo(:,:,:) = (freq(phase,c,n,k)+freq(phase,2,n,k)+freq(phase,3,n,k)+freq(phase,4,n,k)+freq(phase,5,n,k))./nc;
top = ((abs(one(1,1,n,k))).^2.*freq(1,1,n,k))+((abs(one(2,1,n,k))).^2.*freq(2,1,n,k))+(abs(one(3,1,n,k)).^2.*freq(3,1,n,k))+((abs(one(4,1,n,k))).^2.*freq(4,1,n,k));%+(abs((one(5,1,n,k))).^2.*freq(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
bottom =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(one(5,1,n,k))).^2);%+(abs(one(6,1,n,k))).^2);
combo(n,k) = double(top/bottom);

up = ((abs(one(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(one(2,1,n,k))).^2.*pha(2,1,n,k))+(abs(one(3,1,n,k)).^2.*pha(3,1,n,k))+((abs(one(4,1,n,k))).^2.*pha(4,1,n,k));%+(abs((one(5,1,n,k))).^2.*pha(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
down =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(one(5,1,n,k))).^2);%+(abs(one(6,1,n,k))).^2);
together(n,k) = double(up/down);
        end
    end
% phadif = together./(deltate*2*pi);   
figure; map = imagesc(flipud(combo)); colormap(gray);  axis equal; title('GRE FIELD MAP '); set(gca,'Visible','off');hold on
% figure; map = imagesc(flipud(phadif)); colormap(gray);  axis equal; title('GRE FIELD MAP '); set(gca,'Visible','off');hold on
%     plot(x,y, 'r.-', 'MarkerSize', 10)
%     plot(a,b, 'b.-', 'MarkerSize', 10)
%     plot(c,d, 'g.-', 'MarkerSize', 10)
%     
%     LA_filter = double(LA);
% %     section = together.*LA_filter;
% %     for i = 1:npe; for j = 1:ns; if section(i,j) == 0; section(i,j) = finalIm(i,j); end; end; end
% %     imagesc(section)
%     LA_filter(LA_filter == 0) = NaN;
%     leftatrium = combo.*LA_filter;
%     LA_mean = mean(leftatrium(~isnan(leftatrium)));
%     
%     LPV_filter = double(LPV);
%     LPV_filter(LPV_filter == 0) = NaN;
%     leftpv = combo.*LPV_filter;
%     LPV_mean = mean(leftpv(~isnan(leftpv)));
%     
%     RPV_filter = double(RPV);
%     RPV_filter(RPV_filter == 0) = NaN;
%     rightpv = combo.*RPV_filter;
%     RPV_mean = mean(rightpv(~isnan(rightpv)));
%     
%     offset = 0-LA_mean;
%     
%     words = "The mean frequency of the LA is: %0.5f (%d) Hz\nThe mean frequency of the left PV is: %0.5f (%0.5f) Hz\nThe mean frequency of the right PV is: %0.5f (%0.5f) Hz";
%     sprintf(words,LA_mean,0, LPV_mean, LPV_mean +offset, RPV_mean, RPV_mean + offset)
end

% 
% for phase = 1:phases
% figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(one(:,phase,:,:))));
%     colormap(gray); imagesc(flipud((finalIm(:,:)))); axis equal; title('TE = 5ms ');
% 
% figure;finalIm(:,:)=mergeCoilImages(squeeze(abs(two(:,phase,:,:))));
%     colormap(gray); imagesc(flipud((finalIm(:,:)))); axis equal; title('TE = 3.6ms ');hold on 
% 
%     for n = 1:npe
%         for k = 1:ns
% 
% top = ((abs(one(1,1,n,k))).^2.*freq(1,1,n,k))+((abs(one(2,1,n,k))).^2.*freq(2,1,n,k))+(abs(one(3,1,n,k)).^2.*freq(3,1,n,k))+((abs(one(4,1,n,k))).^2.*freq(4,1,n,k));+(abs((one(5,1,n,k))).^2.*freq(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
% bottom =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(one(5,1,n,k))).^2);%+(abs(one(6,1,n,k))).^2);
% combo(n,k) = double(top/bottom);
% 
% up = ((abs(diff(1,1,n,k))).^2.*pha(1,1,n,k))+((abs(diff(2,1,n,k))).^2.*pha(2,1,n,k))+(abs(diff(3,1,n,k)).^2.*pha(3,1,n,k))+((abs(diff(4,1,n,k))).^2.*pha(4,1,n,k));+(abs((diff(5,1,n,k))).^2.*pha(5,1,n,k));%+(abs((one(6,1,n,k))).^2.*freq(6,1,n,k));
% down =  ((abs(one(1,1,n,k))).^2+(abs(one(2,1,n,k))).^2+(abs(one(3,1,n,k))).^2+(abs(one(4,1,n,k))).^2);%+(abs(one(5,1,n,k))).^2);%+(abs(one(6,1,n,k))).^2);
% together(n,k) = double(up/down);
%         end
%     end
% % phadif = together./(deltate*2*pi);   
% figure; map = imagesc(flipud(together)); colormap(gray);  axis equal; title('GRE FIELD MAP '); set(gca,'Visible','off');hold on
% 
% end
