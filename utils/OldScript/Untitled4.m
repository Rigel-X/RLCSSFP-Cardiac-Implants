figure;
for m =1:measure
for j=[1:Nx] 
    for k=[1:Np] 
         if m == 1
             offset = 0;
         end
         if m > 1
             offset = (pi)/(measure*Np);
         end
        mag = (j-Nx/2-1)/Nx; %from -0.5 to 0.5 (without 0.5)
       
        pha=phival(k)+((m-1)*offset);
        loc(m,j,k)=mag*exp(i*pha);

    end
end
plot(real(squeeze(loc(m,:,:))),imag(squeeze(loc(m,:,:))),'c');
hold on 
legend('First Pass: 0 phase increment','Second Pass: 90 phase increment','Third Pass: 180 phase increment','Fourth Pass: 270 phase increment')
end