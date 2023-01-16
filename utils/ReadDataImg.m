clear len a
files=dir('*.IMA');
len=length(files);
for j=1:len
    a(:,:,j)=dicomread(files(j).name);
    figure,imshow(a(:,:,j),[]);%pause(1);
end
img1=double(a(:,:,1));
img2=double(a(:,:,2));
img1 = (img1./max(img1(:))).*2*pi-pi;
img2 = (img2./max(img2(:))).*2*pi-pi;



imfour2=img4_nufft(385:768,384*2+1:384*3,:);
len=size(imfour2,3);
for j = 1:len
imshow(abs(imfour2(:,:,j)),[0,0.0004]);%grid0.04,nufft0.0004,nlinv0.15,irgntv500
title(['phase=',int2str(j),', nufft with DensComp']);colorbar;
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'Sam047_nufft.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'Sam047_nufft.gif','gif','writemode','append','Delaytime',0.05)
end
end




four2=rot90(img_nufft1+img_nufft2+img_nufft3+img_nufft4,2);
len=size(four2,3);
for j = 1:len
    % CV gif
    imshow(abs(four2(:,:,j)),[]);
    title(['phase=',int2str(j)]);colorbar;
    %imshow(abs(four2(97:288,97:288,j)),[]);
    set(gcf,'Color',[1 1 1]);
    set(gcf, 'Position', get(0, 'Screensize'));
    drawnow
    TemGif(j) = getframe(gcf);
    nn=getframe(gcf);
    im=frame2im(nn);
    [I,map]=rgb2ind(im,256);
    if j==1
        imwrite(I,map,'ball068.gif','gif','loopcount',inf,'Delaytime',0.1)
    else
        imwrite(I,map,'ball068.gif','gif','writemode','append','Delaytime',0.1)
    end
end
clear


coilsen=squeeze(four(12,:,:,:));
len=size(coilsen,1);
for j = 1:len
    % CV gif
    imshow(squeeze(abs(coilsen(j,:,:))),[]);title(['phase=',int2str(j)]);axis square,set(gca,'position',[0 0 1 1]),colormap('gray'),colorbar;
    %imshow(four2(:,:,j),[]),title(['phase=',int2str(j)]),axis square,set(gca,'position',[0 0 1 1]),colormap('gray'),colorbar;
    TemGif(j) = getframe(gcf);
    nn=getframe(gcf);
    im=frame2im(nn);
    [I,map]=rgb2ind(im,256);
    if j==1
        imwrite(I,map,'34coils.gif','gif','loopcount',inf,'Delaytime',1)
    else
        imwrite(I,map,'34coils.gif','gif','writemode','append','Delaytime',1)
    end
end


sos4=squeeze(sos(1,:,:,:));
sos4(:,1:384,385:384*2)=squeeze(sos(2,:,:,:));
sos4(:,385:384*2,385:384*2)=squeeze(sos(4,:,:,:));
sos4(:,385:384*2,1:384)=squeeze(sos(3,:,:,:));
sos4(:,385:384*2,384*2+1:384*3)=squeeze(sos(1,:,:,:)+sos(2,:,:,:)+sos(3,:,:,:)+sos(4,:,:,:))/4;






imfour2(1,:,:,:)=finalImg01;
imfour2(2,:,:,:)=finalImg02;
imfour2(3,:,:,:)=finalImg04;
imfour2(4,:,:,:)=finalImg08;
imfour2(5,:,:,:)=finalImg12;
imfour2(6,:,:,:)=finalImg16;
imfour2(7,:,:,:)=finalImg20;
imfour2(8,:,:,:)=finalImg24;
len=size(imfour2,3);
for j = 1:len
if j>=3
    R=4*(j-2);
else
    R=j;
end
imshow(squeeze(abs(imfour2(j,:,:,3))),[]);
title(['acceleration=',int2str(R),', grid']);colorbar;
set(gcf,'Color',[1 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
drawnow
TemGif(j) = getframe(gcf);
nn=getframe(gcf);
im=frame2im(nn);
[I,map]=rgb2ind(im,256);
if j==1
imwrite(I,map,'BL17_grid.gif','gif','loopcount',inf,'Delaytime',0.05)
else
imwrite(I,map,'BL17_grid.gif','gif','writemode','append','Delaytime',0.05)
end
end