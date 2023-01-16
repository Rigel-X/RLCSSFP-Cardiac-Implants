function [recon,cmap]=AdaptiveCombine(im,donorm,rn)
%   Adaptive recon based on Walsh et al.
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery.
%   Magn Reson Med. 2000 May;43(5):682-90.
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob.
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization,
%   Proceedings of the Tenth  Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)
%
%   IN:         im                 image to be reconstructed          (#coils, Ky, Kx)
%               donorm             flag determining whether to normalize image intensity
%               rn                 input noise correlation matrix.
%
%   OUT:        recon              Reconstructed image     ( Ny, Nx)
%               cmap           "Coil maps"           (# coils,  Ny, Nx)
%
%   This non-optimized function will calculate adaptively estimated coil maps
%   based on the Walsh algorithm for use in either optimal combination of
%   array images, or for parallel imaging applications. The donorm flag can be
%   used to include the normalization described in the abstract above. This is
%   only optimal for birdcage type arrays, but will work for many other geometries.
%   The rn matrix should be the noise covariances as decribed in Walsh et al.
%   This could also be a covariance of a region of distrubing artifact as
%   described by Walsh et al.
%
%   The default block size is 4x4. One can also use interpolation to speed
%   up the calculation (see code), but this can cause errors in practice.

%   ~1/1/2001  Mark Griswold
%   FT 2007/03/07

[nc,ny,nx]=size(im);

[mm,maxcoil]=max(sum(sum(permute(abs(im),[3 2 1]))));
maxcoil;



if nargin<3
    rn=eye(nc);
end

if nargin<2
    donorm=0;
end


bs1=4;  %x-block size
bs2=4;  %y-block size
st=1;   %increase to set interpolation step size


wsmall=zeros(nc,round(ny./st),nx./st);
cmapsmall=zeros(nc,round(ny./st),nx./st);


for x=st:st:nx
    for y=st:st:ny
        
        ymin1=max([y-bs1./2 1]);
        xmin1=max([x-bs2./2 1]);
        
        ymax1=min([y+bs1./2 ny]);
        xmax1=min([x+bs2./2 nx]);
        
        ly1=length(ymin1:ymax1);
        lx1=length(xmin1:xmax1);
        
        
        m1=reshape(im(:,ymin1:ymax1,xmin1:xmax1),nc,lx1*ly1);
        
        
        m=m1*m1';
        
        [e,v]=eig(inv(rn)*m);
        v=diag(v);
        [mv,ind]=max(v);
        
        mf=e(:,ind);
        mf=mf/(mf'*inv(rn)*mf);
        normmf=e(:,ind);
        
        mf=mf.*exp(-j*angle(mf(maxcoil)));
        normmf=normmf.*exp(-j*angle(normmf(maxcoil)));
        
        wsmall(:,y./st,x./st)=mf;
        cmapsmall(:,y./st,x./st)=normmf;
        
    end
end

recon=zeros(ny,nx);


for i=1:nc
    wfull(i,:,:)=conj(imresize(squeeze(wsmall(i,:,:)),[ny nx],'bilinear'));
    cmap(i,:,:)=imresize(squeeze(cmapsmall(i,:,:)),[ny nx],'bilinear');
end



wfull=wfull(:,1:ny,1:nx);
cmap=cmap(:,1:ny,1:nx);


for i=1:nc
    recon=recon+squeeze(wfull(i,:,:).*im(i,:,:));
end

if donorm
    recon=recon.*squeeze(sum(abs(cmap).^2));
end
