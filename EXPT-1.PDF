matlab
clear;
clc;
close all;

% Design Parameters
Fs = 8000;
Ap = 1;
As = 40;

% ==========================================
% 1. LOW-PASS FILTER DESIGN
% ==========================================
fp_lp = 1000;
fs_lp = 1500;

Wp_lp = 2*pi*fp_lp/Fs;
Ws_lp = 2*pi*fs_lp/Fs;

ep2 = 10^(Ap/10) - 1;
es2 = 10^(As/10) - 1;

N_lp = ceil(log10(es2/ep2) / ...
    (2*log10(tan(Ws_lp/2)/tan(Wp_lp/2))));

Wn_lp = 2*atan(tan(Wp_lp/2)/(ep2^(1/(2*N_lp))));

% ==========================================
% 2. HIGH-PASS FILTER DESIGN
% ==========================================
fp_hp = 2500;
fs_hp = 1800;

Wp_hp = 2*pi*fp_hp/Fs;
Ws_hp = 2*pi*fs_hp/Fs;

N_hp = ceil(log10(es2/ep2) / ...
    (2*log10(tan(Wp_hp/2)/tan(Ws_hp/2))));

Wn_hp = 2*atan(tan(Wp_hp/2)*(ep2^(1/(2*N_hp))));

% ==========================================
% DESIGN DIGITAL BUTTERWORTH FILTERS
% USING POLE TRANSFORMATION
% ==========================================

% Low-Pass Filter
[b_lp,a_lp] = butterworth_manual(N_lp,Wn_lp,'low');

% High-Pass Filter
[b_hp,a_hp] = butterworth_manual(N_hp,Wn_hp,'high');

% ==========================================
% DISPLAY NUMERICAL OUTPUTS
% ==========================================

fprintf('BUTTERWORTH IIR FILTER DESIGN\n\n');

fprintf('LOW-PASS FILTER\n');
fprintf('Filter Order = %d\n',N_lp);
fprintf('Cutoff Frequency = %.2f Hz\n',Wn_lp*Fs/(2*pi));
fprintf('Numerator Coefficients:\n');
disp(b_lp);
fprintf('Denominator Coefficients:\n');
disp(a_lp);

fprintf('\nHIGH-PASS FILTER\n');
fprintf('Filter Order = %d\n',N_hp);
fprintf('Cutoff Frequency = %.2f Hz\n',Wn_hp*Fs/(2*pi));
fprintf('Numerator Coefficients:\n');
disp(b_hp);
fprintf('Denominator Coefficients:\n');
disp(a_hp);

% ==========================================
% FREQUENCY RESPONSE
% ==========================================

f = linspace(0,Fs/2,1024);
z = exp(1j*2*pi*f/Fs);

H_lp = polyval(fliplr(b_lp),z.^(-1)) ./ ...
       polyval(fliplr(a_lp),z.^(-1));

H_hp = polyval(fliplr(b_hp),z.^(-1)) ./ ...
       polyval(fliplr(a_hp),z.^(-1));

mag_lp = 20*log10(max(abs(H_lp),1e-10));
mag_hp = 20*log10(max(abs(H_hp),1e-10));

phase_lp = angle(H_lp)*180/pi;
phase_hp = angle(H_hp)*180/pi;

% ==========================================
% PLOTS
% ==========================================

figure('Name','Butterworth IIR Filters');

subplot(2,2,1);
plot(f,mag_lp,'b','LineWidth',1.5);
grid on;
title('Low-Pass Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-80 5]);

subplot(2,2,2);
plot(f,phase_lp,'r','LineWidth',1.5);
grid on;
title('Low-Pass Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

subplot(2,2,3);
plot(f,mag_hp,'b','LineWidth',1.5);
grid on;
title('High-Pass Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-80 5]);

subplot(2,2,4);
plot(f,phase_hp,'r','LineWidth',1.5);
grid on;
title('High-Pass Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

% ==========================================
% LOCAL FUNCTION
% ==========================================

function [b,a] = butterworth_manual(N,Wn,type)

    % Prewarp the digital cutoff frequency
    Wc = 2*tan(Wn/2);

    % Analog Butterworth poles
    k = 1:N;
    poles = Wc*exp(1j*pi*(2*k+N-1)/(2*N));

    % Bilinear transformation
    zp = (2+poles)./(2-poles);

    % Digital filter coefficients
    a = real(poly(zp));

    if strcmp(type,'low')
        % Low-pass zeros at z = -1
        zzeros = -ones(1,N);
    else
        % High-pass zeros at z = 1
        zzeros = ones(1,N);
    end

    b = real(poly(zzeros));

    % Normalize the gain
    if strcmp(type,'low')
        gain = sum(a)/sum(b);
    else
        gain = sum(a.*((-1).^(0:N)))/ ...
               sum(b.*((-1).^(0:N)));
    end

    b = b*gain;
end

