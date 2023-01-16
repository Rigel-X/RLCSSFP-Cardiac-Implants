%true ram-lak filter, Thomas Thuring, 30/04/07

function coeff = trueRamLak(N)

if(mod(N/2,2) == 0) %check if N/2 is even
    Nr=N+2;
else
    Nr=N;
end

%create ranges
n=[1:Nr];
odd=[1:2:Nr];
odd2=[-Nr/2:2:-1,1:2:Nr/2-2];

%create coefficients in spatial domain
f=zeros(1,Nr);
f(odd)=-2./(pi*(odd2*1).^2); %dt=1
f(Nr/2+1)=pi/(2*(1)^2);
plot(f);

%fft to get k-space coeffs
k=abs(fftshift(fft(f)));
k2=k./(2*max(k));
coeff=k([2:N+1]);
