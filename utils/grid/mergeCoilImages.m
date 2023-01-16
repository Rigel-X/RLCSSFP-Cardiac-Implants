% merge images obtained from different coils, Thomas Thuring, 27/04/07
% coildata: 3D-array of form imageData(coilnumber, Ny, Nx)
% with Ny: # of rows in image
%      Nx: # of columns image
% function takes data of all coils at one location and computes output
% image as follows: out=sqrt(sum(coildatapoint^2,for all coils))

function image = mergeCoilImages(imageData)

if(size(size(imageData),2)==2)  
    data(1,:,:)=imageData; %algorithm takes only 3 dimensional
else
    data=imageData;
end

Ncoils = size(data,1);
Ny = size(data,2);
Nx = size(data,3);

temp=zeros(Ny,Nx);
for coil=1:Ncoils  % for all coil elements
    coilImage = squeeze(data(coil,:,:));
    temp = temp+(coilImage.^2);
end
image = sqrt(temp);