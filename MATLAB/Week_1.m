Fs = 8000; % Sampling rate
T = 0.2;   % Duration

Fm = 10;   % Message signal frequency
Ac = 1;    % Carrier amplitude
Fc = 300;  % Carrier frequency

samples = Fs * T;

% Time vector
t = linspace(0, T, samples + 1); t(end) = [];
% Frequency vector
f = linspace(-Fs / 2, Fs / 2, samples + 1); f(end) = [];

x = cos(2 * pi * Fm * t);      % Message signal
c = Ac * cos(2 * pi * Fc * t); % Carrier waveform

%% DSB-SC modulation
y1 = x .* c;

figure(), title('DSB-SC'), subplot(2, 1, 1)
plot(t, x, 'm', 'LineWidth', 1), hold on
plot(t, c, 'c--'), plot(t, y1, 'k--')
xlabel('Time (sec)'), ylabel('Amplitude')
legend('m(t)', 'c(t)', 's(t)')

% Spectral Characteristics of the DSB-SC signal
X = fft(x); C = fft(c); Y = fft(y1);

subplot(2, 1, 2)
stem(f, abs(fftshift(X)) / Fs, '.m', 'LineWidth', 1), hold on
stem(f, abs(fftshift(C)) / Fs, '.c')
stem(f, abs(fftshift(Y)) / Fs, '.k')
xlabel('Frequency (Hz)'), ylabel('Magnitude')
legend('M(f)', 'C(f)', 'Y(f)'), axis([-2 * Fc 2 * Fc 0 0.15])

%% DSB-FC modulation
mu = 2; % Modulating index
y2 = (1 + mu * x) .* c; 

figure(), title('DSB-FC'), subplot(2, 1, 1)
plot(t, x, 'm', 'LineWidth', 1), hold on
plot(t, c, 'c--'), plot(t, y2, 'k--')
xlabel('Time (sec)'), ylabel('Amplitude')
legend('m(t)', 'c(t)', 'y(t)')

% Task - Change the modulation index, and observe the resultant modulated signal.
%        Comment on the results

% Spectral Characteristics of the DSB-FC signal
Y = fft(y2);

subplot(2, 1, 2)
stem(f, abs(fftshift(X)) / Fs, '.m', 'LineWidth', 1), hold on
stem(f, abs(fftshift(C)) / Fs, '.c')
stem(f, abs(fftshift(Y)) / Fs, '.k')
xlabel('Frequency (Hz)'), ylabel('Magnitude')
legend('M(f)', 'C(f)', 'Y(f)'), axis([-2 * Fc 2 * Fc 0 0.15])

%% Recovering the original signal
x_baseband = y2;
x_baseband(x_baseband < 0) = 0; % Rectified Signal

X_filt = fft(x_baseband);               % Spectral components
X_filt(abs(ifftshift(f)) > 2 * Fm) = 0; % Filtering to remove high frequency carrier

x_demod = real(ifft(X_filt));      % Recover the signal
x_demod = x_demod - mean(x_demod); % Remove the DC bias

% Plot against original
figure()
hold on
plot(t, x, 'k')
plot(t, x_demod, 'r--')
xlabel('Time [s]'), ylabel('Amplitude')
legend('x(t)', 'x_{rx}(t)')
