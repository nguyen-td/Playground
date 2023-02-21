% Compute discrete time Fourier transform and compare it to fast Fourier
% transform.

%% Create a mixed signal
freq1 = 3; % frequency in Hz
freq2 = 8; % frequency in Hz

fs = 100;        % samples per second
dt = 1/fs;       % seconds per sample
T = 1;           % signal length in seconds
t = (0:dt:T-dt); % array of time points
N = fs*T;        % number of samples/data points

[mixed, sinusoids] = create_sinusoid(fs,T,{freq1, freq2});
sin1 = sinusoids{1};
sin2 = sinusoids{2};

%% Plot sine waves
subplot(2,3,1)
plot(sin1);
hold on;
plot(sin2);
xlabel("time")
title("3 and 8 Hz oscillations")
legend("3 Hz signal", "8 Hz signal")

subplot(2,3,4)
plot(mixed)
xlabel("time")
title('Mixed signal')

%% Run manual DTFT
y_dtft = zeros(1, fs); % initialize Fourier coefficients
for fi=1:fs
    % create sine wave
    sine_wave = exp(-1i*2*pi*(fi-1).*t); % . creates element-wise operations
    % compute dot product between sine wave and mixed signal
    y_dtft(fi) = sum(sine_wave.*mixed);
end
y_dtfn = y_dtft/N;
%% Run FFT
y_fft = fft(mixed);
% f = (0:length(y_fft)-1)*Fs/length(y_fft);
%% Plot frequency components
subplot(2,3,2)
plot(abs(y_fft).^2)
title("Power spectrum (FFT)")
xlabel("Frequency [Hz]")
ylabel("Power")

subplot(2,3,5)
plot(abs(y_dtft).^2)
title("Power spectrum (DTFT)")
xlabel('Frequency [Hz]')
ylabel('Power')

%% Keep only positivive frequencies
nyquist = N/2+1;
y_fft_ny = y_fft(1:nyquist);
y_dtft_ny = y_dtft(1:nyquist);

subplot(2,3,3)
plot(abs(y_fft_ny).^2)
title("PSD without negative frequencies (FFT)")
xlabel("Frequency [Hz]")
ylabel("Power")

subplot(2,3,6)
plot(abs(y_dtft_ny).^2)
title("PSD without negative frequencies (DTFT)")
xlabel('Frequency [Hz]')
ylabel('Power')

%% Sampling below Nyquist (try example in Steve Brunton)
