
% for m = 1:1
    
%          if m == 1
%              offset = 0;
%          end
%          if m > 1
             offset = (pi)/(measure*Np);
%          end
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
       
        pha=phival(k)+((m-1)*offset);
        loc(j,k)=mag*exp(i*pha);
        
% C = ['c' 'b' 'g' 'm'];
location = squeeze(loc(:,:));
h = plot(real(location),imag(location), 'color', 'b');
hc = get(h, 'Children');
% hold on 
% h(m) = plot(real(squeeze(loc(m,:,:))),imag(squeeze(loc(m,:,:))), 'color', C(m));hold on 
legend(['First Pass: 0 phase increment','Second Pass: 90 phase increment','Third Pass: 180 phase increment','Fourth Pass: 270 phase increment')
