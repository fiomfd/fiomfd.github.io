% computer_sound_a.m
%%%%%%%%%%%%%%%%%%%%%%
%
% Create a signal
%
%%%%%%%%%%%%%%%%%%%%%%
fs=1000;
Ts=1/fs;
N=12000;
t=(0:N-1)*Ts;
s=sin(2*t.*(t-3).*(t-6).*(t-9).*(t-12));
% sound(s)
% plot(t,s);
% title('Original Audio Signal x(t)');
% xlabel('time t [sec]');
% ylabel('amplitude');

plot(s);
title('Discrete audio signal x[n]:=x(n/fs)','x(t) = sin(2t(t-3)(t-6)(t-9)(t-12))');
xlabel('n = 0, 1, 2, ..., N-1');
ylabel('x[n]');
yticks([-1 0 1]);
