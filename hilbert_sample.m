fs = 10000;
t = 0:1/fs:2-1/fs;
x = sin(2*pi*60*t);

z = hilbert(x);
y = x + i*z;

figure
plot(t2,real(y2),t2,imag(y2))
xlim([0.01 0.03])              
legend('real','imaginary')     
title('FIR Filter')



instfrq = 1000/(2*pi)*diff(unwrap(angle(z)));

figure
plot(t(2:end),instfrq)
ylim([0 10])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
