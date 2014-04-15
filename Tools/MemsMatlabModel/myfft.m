function [f, amp] = myfft(data, tsamp)

Fs = 1/tsamp;
T = tsamp;
L = length(data);

x = data;

NFFT = 2^nextpow2(L);
Y = fft(x,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

amp = 2*abs(Y(1:NFFT/2+1));
loglog(f,amp)
title('Single-Sided Amplitude Spectrum of data(t)')
xlabel('Frequency (Hz)')
ylabel('|Data(f)|')
end