% computer_sound.m
%%%%%%%%%%%%%%%%%%%%%%
%
% Combine all parts
% - Create signal
% - discrete Fourier transform
% - spectrogram
%%%%%%%%%%%%%%%%%%%%%%
fs=1000;
Ts=1/fs;
L=12000;
t=(0:L-1)*Ts;
s=sin(2*t.*(t-3).*(t-6).*(t-9).*(t-12));

hat = fft(s);

wft=spectrogram(s);

subplot(3,1,1);
plot(t,s);
title('Original Audio Signal x(t)');
xlabel('time t [sec]');
ylabel('amplitude');

subplot(3,1,2); 
fshift = (-L/2:L/2-1)*(100/L);
hatshift = fftshift(hat);
plot(fshift,abs(hatshift));
title('Discrete Fourier transform');
xlabel('discrete frequency k');
ylabel('|F(x)[k]|');

subplot(3,1,3);
spectrogram(s,'yaxis');
title('Spectrogram of Audio Signal');
