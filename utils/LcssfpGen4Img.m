% jie xiang @yale mrrc
% generate the four bssfp images of lcssfp to check different bands pattern

disp('------------------- image1 -------------------------');
if isphaseinterleaved == 1
for ncoil = 1:coils
    for xi = 1:Nx
        for yi = 1:Nx
            TimeCurve(:,ncoil,xi,yi) = imone(:,ncoil,xi,yi);
            % TimeCurve = abs(four(:,ncoil,xi,yi));
            FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,ncoil,xi,yi))))).*(CropWindow');
            imone(:,ncoil,xi,yi)=fftshift(ifft(ifftshift(FreCurve)));
            %four(:,ncoil,xi,yi)=FreCurve;
        end
    end
end
end
imone2=zeros(Nx,Nx,phases);
% for ph = 1:phases
%     imone2(:,:,ph) = AdaptiveCombine(squeeze(imone(ph,:,:,:)));
% end
imone2 = squeeze(sqrt(sum(abs(imone).^2,2)));
imone2 = permute(imone2,[2,3,1]);

% len=size(imone2,3);
% for j = 1:len
% % CV gif
% %imshow(abs(imone2(:,:,j)),[]);
% imshow(abs(imone2(:,:,j)),[0,0.05]);
% title(['phase=',int2str(j),', slice=4']);colorbar;
% %title(['phase=',int2str(j)]);colorbar;
% %imshow(abs(four2(97:288,97:288,j)),[]);
% set(gcf,'Color',[1 1 1]);
% set(gcf, 'Position', get(0, 'Screensize'));
% drawnow
% TemGif(j) = getframe(gcf);
% nn=getframe(gcf);
% im=frame2im(nn);
% [I,map]=rgb2ind(im,256);
% if j==1
% imwrite(I,map,'Sam048_1.gif','gif','loopcount',inf,'Delaytime',0.05)
% else
% imwrite(I,map,'Sam048_1.gif','gif','writemode','append','Delaytime',0.05)
% end
% end
% close all

disp('------------------- image2 -------------------------');
if isphaseinterleaved == 1
for ncoil = 1:coils
    for xi = 1:Nx
        for yi = 1:Nx
            TimeCurve(:,ncoil,xi,yi) = imtwo(:,ncoil,xi,yi);
            % TimeCurve = abs(four(:,ncoil,xi,yi));
            FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,ncoil,xi,yi))))).*(CropWindow');
            imtwo(:,ncoil,xi,yi)=fftshift(ifft(ifftshift(FreCurve)));
            %four(:,ncoil,xi,yi)=FreCurve;
        end
    end
end
end
imtwo2=zeros(Nx,Nx,phases);
% for ph = 1:phases
%     imtwo2(:,:,ph) = AdaptiveCombine(squeeze(imtwo(ph,:,:,:)));
% end
imtwo2 = squeeze(sqrt(sum(abs(imtwo).^2,2)));
imtwo2 = permute(imtwo2,[2,3,1]);
% 
% len=size(imtwo2,3);
% for j = 1:len
% % CV gif
% %imshow(abs(imone2(:,:,j)),[]);
% imshow(abs(imtwo2(:,:,j)),[0,0.05]);
% title(['phase=',int2str(j),', slice=4']);colorbar;
% %title(['phase=',int2str(j)]);colorbar;
% %imshow(abs(four2(97:288,97:288,j)),[]);
% set(gcf,'Color',[1 1 1]);
% set(gcf, 'Position', get(0, 'Screensize'));
% drawnow
% TemGif(j) = getframe(gcf);
% nn=getframe(gcf);
% im=frame2im(nn);
% [I,map]=rgb2ind(im,256);
% if j==1
% imwrite(I,map,'Sam048_2.gif','gif','loopcount',inf,'Delaytime',0.05)
% else
% imwrite(I,map,'Sam048_2.gif','gif','writemode','append','Delaytime',0.05)
% end
% end
% close all


disp('------------------- image3 -------------------------');
if isphaseinterleaved == 1
for ncoil = 1:coils
    for xi = 1:Nx
        for yi = 1:Nx
            TimeCurve(:,ncoil,xi,yi) = imthree(:,ncoil,xi,yi);
            % TimeCurve = abs(four(:,ncoil,xi,yi));
            FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,ncoil,xi,yi))))).*(CropWindow');
            imthree(:,ncoil,xi,yi)=fftshift(ifft(ifftshift(FreCurve)));
            %four(:,ncoil,xi,yi)=FreCurve;
        end
    end
end
end
imthree2=zeros(Nx,Nx,phases);
% for ph = 1:phases
%     imthree2(:,:,ph) = AdaptiveCombine(squeeze(imthree(ph,:,:,:)));
% end
imthree2 = squeeze(sqrt(sum(abs(imthree).^2,2)));
imthree2 = permute(imthree2,[2,3,1]);

% len=size(imthree2,3);
% for j = 1:len
% % CV gif
% %imshow(abs(imone2(:,:,j)),[]);
% imshow(abs(imthree2(:,:,j)),[0,0.05]);
% title(['phase=',int2str(j),', slice=4']);colorbar;
% %title(['phase=',int2str(j)]);colorbar;
% %imshow(abs(four2(97:288,97:288,j)),[]);
% set(gcf,'Color',[1 1 1]);
% set(gcf, 'Position', get(0, 'Screensize'));
% drawnow
% TemGif(j) = getframe(gcf);
% nn=getframe(gcf);
% im=frame2im(nn);
% [I,map]=rgb2ind(im,256);
% if j==1
% imwrite(I,map,'Sam048_3.gif','gif','loopcount',inf,'Delaytime',0.05)
% else
% imwrite(I,map,'Sam048_3.gif','gif','writemode','append','Delaytime',0.05)
% end
% end
% close all

disp('------------------- image4 -------------------------');
if isphaseinterleaved == 1
for ncoil = 1:coils
    for xi = 1:Nx
        for yi = 1:Nx
            TimeCurve(:,ncoil,xi,yi) = imfour(:,ncoil,xi,yi);
            % TimeCurve = abs(four(:,ncoil,xi,yi));
            FreCurve = ifftshift(fft(fftshift(squeeze(TimeCurve(:,ncoil,xi,yi))))).*(CropWindow');
            imfour(:,ncoil,xi,yi)=fftshift(ifft(ifftshift(FreCurve)));
            %four(:,ncoil,xi,yi)=FreCurve;
        end
    end
end
end
imfour2=zeros(Nx,Nx,phases);
% for ph = 1:phases
%     imfour2(:,:,ph) = AdaptiveCombine(squeeze(imfour(ph,:,:,:)));
% end
imfour2 = squeeze(sqrt(sum(abs(imfour).^2,2)));
imfour2 = permute(imfour2,[2,3,1]);
% 
% len=size(imfour2,3);
% for j = 1:len
% % CV gif
% %imshow(abs(imone2(:,:,j)),[]);
% imshow(abs(imfour2(:,:,j)),[0,0.05]);
% title(['phase=',int2str(j),', slice=4']);colorbar;
% %title(['phase=',int2str(j)]);colorbar;
% %imshow(abs(four2(97:288,97:288,j)),[]);
% set(gcf,'Color',[1 1 1]);
% set(gcf, 'Position', get(0, 'Screensize'));
% drawnow
% TemGif(j) = getframe(gcf);
% nn=getframe(gcf);
% im=frame2im(nn);
% [I,map]=rgb2ind(im,256);
% if j==1
% imwrite(I,map,'Sam048_4.gif','gif','loopcount',inf,'Delaytime',0.05)
% else
% imwrite(I,map,'Sam048_4.gif','gif','writemode','append','Delaytime',0.05)
% end
% end

% disp('------------------- image5 -------------------------');
% for j = 1:len
% subplot(2,2,1)
% imshow(rot90(abs(imone2(:,:,j))),[0,0.05]);
% title(['phase=',int2str(j),', slice=1, phase cycle 1']);
% subplot(2,2,2)
% imshow(rot90(abs(imtwo2(:,:,j))),[0,0.05]);
% title(['phase=',int2str(j),', slice=1, phase cycle 2']);
% subplot(2,2,3)
% imshow(rot90(abs(imthree2(:,:,j))),[0,0.05]);
% title(['phase=',int2str(j),', slice=1, phase cycle 3']);
% subplot(2,2,4)
% imshow(rot90(abs(imfour2(:,:,j))),[0,0.05]);
% title(['phase=',int2str(j),', slice=1, phase cycle 4']);
% set(gcf,'Color',[1 1 1]);
% set(gcf, 'Position', get(0, 'Screensize'));
% drawnow
% TemGif(j) = getframe(gcf);
% nn=getframe(gcf);
% im=frame2im(nn);
% [I,map]=rgb2ind(im,256);
% if j==1
% imwrite(I,map,'Sam048_.gif','gif','loopcount',inf,'Delaytime',0.05)
% else
% imwrite(I,map,'Sam048_.gif','gif','writemode','append','Delaytime',0.05)
% end
% end
img=zeros(Nx,5*Nx,phases);
img(:,1:Nx,:)=imone2(:,:,:);
img(:,Nx+1:2*Nx,:)=imtwo2(:,:,:);
img(:,2*Nx+1:3*Nx,:)=imthree2(:,:,:);
img(:,3*Nx+1:4*Nx,:)=imfour2(:,:,:);
img(:,4*Nx+1:5*Nx,:)=four2(:,:,:);%./3;
figure,imshow(imresize(squeeze(abs(img(:,:,3))),2),[]);

%save('Recon_Image\Dana\DP39_4d_pd60_20p_acc_ave.mat','imone2','imtwo2','imthree2','imfour2','img','-append');


if isdynamicinterleaved==0
lcSSFP1 = imone2;lcSSFP2 = imtwo2;lcSSFP3 = imthree2;lcSSFP4 = imfour2;lcSSFP0 = img;
save('ksptest.mat','lcSSFP1','lcSSFP2','lcSSFP3','lcSSFP4','lcSSFP0','-append');
end

if isdynamicinterleaved==1 && isphaseinterleaved == 0
idSSFP1 = imone2;idSSFP2 = imtwo2;idSSFP3 = imthree2;idSSFP4 = imfour2;idSSFP0 = img;
save('ksptest.mat','idSSFP1','idSSFP2','idSSFP3','idSSFP4','idSSFP0','-append');
end


if isdynamicinterleaved==1 && isphaseinterleaved == 1
ipSSFP1 = imone2;ipSSFP2 = imtwo2;ipSSFP3 = imthree2;ipSSFP4 = imfour2;ipSSFP0 = img;
save('ksptest.mat','ipSSFP1','ipSSFP2','ipSSFP3','ipSSFP4','ipSSFP0','-append');
end




