% Compare reconstruction by complying to the Nyquist-Shannon sampling 
% theorem, i.e., by sampling at least at twice the highest frequency, to
% reconstruction by violating against the theorem, i.e., sampling below the
% Nyquist frequency.

%% Create a mixed signal at Nyquist frequency
freq1 = 20; % frequency in Hz
freq2 = 40; % frequency in Hz
nyquist = 2*freq2+1;

fs = nyquist;   % samples per second
dt = 1/fs;       % seconds per sample
T = 1;           % signal length in seconds
t = (0:dt:T-dt); % array of time points
N = fs*T;        % number of samples/data points

[mixed_atnyq, sinusoids] = create_sinusoid(fs,T,{freq1, freq2});
sin1 = sinusoids{1};
sin2 = sinusoids{2};

%% Plot sine waves
subplot(2,3,1)
plot(sin1);
hold on;
plot(sin2);
xlabel("time")
title("Sampled at Nyquist frequency")
legend("20 Hz", "40 Hz")

subplot(2,3,4)
plot(mixed_atnyq)
xlabel("time")

%% Create a mixed signal below Nyquist frequency
fs = nyquist-10;   % samples per second
dt = 1/fs;       % seconds per sample
T = 1;           % signal length in seconds

[mixed_beloqnyq, sinusoids] = create_sinusoid(fs,T,{freq1, freq2});
sin1 = sinusoids{1};
sin2 = sinusoids{2};

%% Plot sine waves
subplot(2,3,2)
plot(sin1);
hold on;
plot(sin2);
xlabel("time")
title("Sampled below Nyquist frequency")
legend("20 Hz", "40 Hz")

subplot(2,3,5)
plot(mixed_beloqnyq)
xlabel("time")

%% Run FFT
y_fft_atnyq = fft(mixed_atnyq);
y_fft_beloqnyq = fft(mixed_beloqnyq);

subplot(2,3,3)
plot(abs(y_fft_atnyq).^2)
title("PSD (signal at Nyquist frequency)")
xlabel("Frequency [Hz]")
ylabel("Power")

subplot(2,3,6)
plot(abs(y_fft_beloqnyq).^2)
title("PSD (signal below Nyquist frequency)")
xlabel('Frequency [Hz]')
ylabel('Power')