

dir1='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1019_DP\timtriob_20181019_071223932\DP.S011.D500181019\'
dir3='C:\Users\Ricardo\Desktop\JaimesFiles\Volunteers\1019_DP\timtriob_20181019_071223932\DP.S012.D500181019\'
[header , volumecar]=dicomvolumeget(dir1);
[header , volumerad]=dicomvolumeget(dir3)

lccine=squeeze(final(:,1,:,:)); %then normalize each.

carcine=volumecar(:,:,1:1:end);
carcine2(21:140,1:160,1:20)=carcine;
radcine=volumerad(:,:,1:1:end);
figure(77);
hold on;caxis([0 1000]);
%%
v = VideoWriter('lcssfpmoviesnew.avi');

v.FrameRate = 5;
 open(v)

clear im3;
minval=min(min(min(carcine(:,:,:))));
maxval=max(max(max(carcine(:,:,:))));

maxvalr=max(max(max(radcine(:,:,:))));
maxvall=max(max(max(lccine(:,:,:))));

radcine1=radcine*(maxval/maxvalr);

lccine1=lccine*(maxval/maxvall);
for i=1:15,
 %not flipupd because not correct now wiht imshow.
    im1=squeeze((carcine2(:,:,i)));
 
     im2=squeeze((radcine1(:,:,i)));
    
     im3=squeeze((lccine1(i,:,:)));
    
        composite=[im1,im2,im3];
        size(composite)
        
%            newim1= insertText(composite,[80,10],'Cartesian','TextColor','white', 'BoxColor','black');
%              newim2=insertText(newim1,[240,10],'Radial','TextColor','white', 'BoxColor','black');
%                   newim3=insertText(newim2,[400,10],'MA-SSFP','TextColor','white', 'BoxColor','black');
        imshow(composite,[minval maxval*.5]);% axis equal;axis off; colormap gray; title(dirname);
%                insertText(composite,[80,10],'Cartesian','TextColor','white');
%               insertText(composite,[240,10],'Radial','TextColor','white');
%                  insertText(composite,[400,10],'MA-SSFP','TextColor','white');
   frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(v, frame);

end
close(v)
