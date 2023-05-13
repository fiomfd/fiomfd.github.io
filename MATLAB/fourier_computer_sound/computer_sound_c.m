% computer_sound_c.m
%%%%%%%%%%%%%%%%%%%%%%
fs=1000;
Ts=1/fs;
L=12000;
t=(0:L-1)*Ts;
s=sin(2*t.*(t-3).*(t-6).*(t-9).*(t-12));
%%%%%%%%%%%%%%%%%%%%
%
% spectrogram 
%
%%%%%%%%%%%%%%%%%%%%
wft=spectrogram(s);
spectrogram(s,'yaxis');
title('Spectrogram of Audio Signal');
