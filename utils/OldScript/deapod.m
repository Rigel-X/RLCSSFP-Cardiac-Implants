function [ im ] = deapod( orig,oversmpl,res, kwidth)
%DEAPODIZATION - function computed and applies deapodisation to 2D image array
% im = output image, deapodized image
% orig = original image matrix
% res = resolution of image, size of image matrix i.e. 160, 256
% kwidth = kernel width


osf=oversmpl;
n=res;
kw=kwidth;
beta= pi*sqrt((kw/osf*(osf-0.5)).^2-0.8);

% compute deappodization function
x = [-n/2:n/2-1]/n;
% size(x)
% size(orig)
% inverse fourier transform of kaiser bessel function
% sqa = sqrt(pi*pi*kw*kw*x.*x-beta*beta);
% dax = sin(sqa)./(sqa);
% % normalize by DC value
% dax = dax/dax(n/2+1);
% make it a 2D array

[x1,x2] = meshgrid(x);
dist = sqrt(x1.^2+x2.^2);
dax = sqrt(pi*pi*kw*kw*dist.*dist - beta^2);
dax = sin(dax)./dax;
da = dax/max(dax(:));

weight = dist>(0.5);
% deappodize
im = orig./da;
im(weight) = 0;
end

