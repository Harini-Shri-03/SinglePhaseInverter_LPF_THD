clc; clear; close all;

%% ===== PARAMETERS =====
Vdc = 100;           % DC bus voltage (V)
f_out = 50;          % desired output frequency (Hz)
fs = 10000;          % sampling frequency (Hz) - increase for more accuracy
numCycles = 20;      % integer number of output cycles (helps FFT resolution)
Tsim = numCycles / f_out;
t = 0:1/fs:Tsim-1/fs;    % time vector
N = length(t);

% SPWM parameters
f_carrier = 1000;    % carrier (triangular) frequency Hz
ma = 0.9;            % modulation index (0..1)

%% ===== Generate reference sine and triangular carrier =====
Vref = ma * sin(2*pi*f_out*t);          % reference sine [-ma,ma]

% Triangular carrier (amplitude +/-1). Uses mod() to avoid toolbox dependency.
tri = 4*abs(mod(f_carrier * t, 1) - 0.5) - 1;   % range [-1,1]

%% ===== Generate PWM and Inverter Output (±Vdc) =====
PWM = double(Vref >= tri);         % comparator -> 0/1
Vout = Vdc * (2*PWM - 1);          % map to ±Vdc

%% ===== Plot first few cycles (time-domain) =====
plotCycles = 3;
idxEnd = find(t >= plotCycles/f_out, 1) - 1;

figure('Name','Time-domain signals','NumberTitle','off');
subplot(3,1,1);
plot(t(1:idxEnd), Vref(1:idxEnd), 'LineWidth', 1.1);
title('Reference sine (first few cycles)'); xlabel('Time (s)'); ylabel('V_{ref}');

subplot(3,1,2);
plot(t(1:idxEnd), tri(1:idxEnd), 'LineWidth', 1.1);
title('Triangular carrier (first few cycles)'); xlabel('Time (s)'); ylabel('Carrier');

subplot(3,1,3);
stairs(t(1:idxEnd), Vout(1:idxEnd), 'LineWidth', 1);
title('Inverter PWM output (first few cycles)'); xlabel('Time (s)'); ylabel('V_{out}');

%% ===== FFT / Spectrum and THD (before filter) =====
% Windowing to reduce leakage
win = ones(1,N);         % row vector
v = Vout .* win;

% Zero-pad for better frequency resolution
NFFT = 2^nextpow2(N) * 4;
Y = fft(v, NFFT);

% single-sided amplitude spectrum
P2 = abs(Y)/NFFT;
P1 = P2(1:NFFT/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f1 = (0:(NFFT/2)) * (fs / NFFT);

% find fundamental index (closest bin)
[~, idxFund] = min(abs(f1 - f_out));
A1 = P1(idxFund);

% harmonic indices and THD calculation
maxHarm = floor((fs/2)/f_out);                % up to Nyquist
harmNums = 2:maxHarm;
idxs = round(harmNums * f_out * NFFT / fs) + 1; % +1 for MATLAB indexing
idxs = idxs(idxs <= length(P1));              % discard out-of-range

A_harm = P1(idxs);
THD_before = 100 * sqrt(sum(A_harm.^2)) / A1;  % percent

fprintf('THD (before filter) = %.3f %%\n', THD_before);

% Plot spectrum (before filter)
figure('Name','Spectrum before filter','NumberTitle','off');
plot(f1, P1, 'LineWidth', 1);
xlim([0, min(5000, fs/2)]);
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Single-sided Spectrum - before filter');

%% ===== Design LC low-pass filter =====
Rload = 10;          % load resistance in ohms
fc_cut = 200;        % desired cutoff (Hz) -- set > f_out but << f_carrier
L = 5e-3;            % choose L (H) - example 5 mH
C = 1 / ((2*pi*fc_cut)^2 * L);  % compute C to get approx cutoff

fprintf('LC filter chosen: L = %.3e H, C = %.3e F, fc (design) = %.1f Hz\n', L, C, fc_cut);

%% ===== Apply LC filter in frequency domain =====
% Work with shifted FFT so we can evaluate H(f) from negative to positive freqs
Yshift = fftshift(Y);
fshift = (-NFFT/2 : NFFT/2-1) * (fs / NFFT);    % frequency axis for shifted FFT
s = 1j*2*pi*fshift;

% impedances
Zc = 1 ./ (s * C);          % capacitor impedance (inf at DC)
Zl = s * L;                 % inductor impedance
Zload = 1 ./ (1/Rload + 1 ./ Zc);   % parallel C||R

% transfer function Vout/Vin (series L then parallel C||R)
Hshift = Zload ./ (Zl + Zload);

% apply filter in frequency domain and inverse FFT
Yfilt_shift = Yshift .* Hshift;
Yfilt = ifftshift(Yfilt_shift);
v_filt = ifft(Yfilt, 'symmetric');   % filtered time-domain signal (real)

%% ===== THD after filtering =====
vwin = v_filt(1:N) .* win;
Y2 = fft(vwin, NFFT);
P2f = abs(Y2)/NFFT;
P1f = P2f(1:NFFT/2+1);
P1f(2:end-1) = 2*P1f(2:end-1);

A1f = P1f(idxFund);
A_harm_f = P1f(idxs);
THD_after = 100 * sqrt(sum(A_harm_f.^2)) / A1f;

fprintf('THD (after LC filter) = %.3f %%\n', THD_after);

% Plot spectrum after filter
figure('Name','Spectrum after LC filter','NumberTitle','off');
plot(f1, P1f, 'LineWidth', 1);
xlim([0, min(2000, fs/2)]);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Single-sided Spectrum - after LC filter');

% Plot time domain comparison (zoomed)
figure('Name','Time-domain comparison','NumberTitle','off');
subplot(2,1,1);
stairs(t(1:idxEnd), Vout(1:idxEnd), 'LineWidth', 1);
title('PWM output (time domain, unfiltered)'); xlabel('Time (s)'); ylabel('V');

subplot(2,1,2);
plot(t(1:idxEnd), v_filt(1:idxEnd), 'LineWidth', 1.1);
title(sprintf('Filtered output (time domain) — THD: %.3f %%', THD_after));
xlabel('Time (s)'); ylabel('V');

%% ===== Optional: Sweep modulation index and show THD change =====
mas = 0.2:0.1:1.0;
THD_b = zeros(size(mas));
THD_a = zeros(size(mas));

for k=1:length(mas)
    Vref_k = mas(k) * sin(2*pi*f_out*t);
    PWM_k = double(Vref_k >= tri);
    Vout_k = Vdc * (2*PWM_k - 1);
    vk = Vout_k .* win;
    Yk = fft(vk, NFFT);
    P2k = abs(Yk)/NFFT; P1k = P2k(1:NFFT/2+1); P1k(2:end-1)=2*P1k(2:end-1);
    A1k = P1k(idxFund);
    idxs_k = idxs(idxs <= length(P1k));
    THD_b(k) = 100 * sqrt(sum(P1k(idxs_k).^2)) / A1k;

    % filter
    Yk_shift = fftshift(Yk);
    Yf_k_shift = Yk_shift .* Hshift;
    Yf_k = ifftshift(Yf_k_shift);
    v_filt_k = ifft(Yf_k, 'symmetric');
    v_filt_k_win = v_filt_k(1:N) .* win;
    Yfk = fft(v_filt_k_win, NFFT);
    P2fk = abs(Yfk)/NFFT; P1fk = P2fk(1:NFFT/2+1); P1fk(2:end-1)=2*P1fk(2:end-1);
    A1fk = P1fk(idxFund);
    THD_a(k) = 100 * sqrt(sum(P1fk(idxs_k).^2)) / A1fk;
end

figure('Name','THD vs modulation index','NumberTitle','off');
plot(mas, THD_b, '-o', mas, THD_a, '-x', 'LineWidth', 1.1);
xlabel('Modulation index (m_a)'); ylabel('THD (%)');
legend('Before filter','After filter'); grid on;
title('THD variation with modulation index');

% End of script
