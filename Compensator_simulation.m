clc;

%% ---------------importing data---------------------------

G_freq = readmatrix('Noise_cancelling Headphone data.xlsx','Sheet','Data','Range','A2:A40');
G_mag_dB = readmatrix('Noise_cancelling Headphone data.xlsx','Sheet','Data','Range','L2:L40');
G_phase = readmatrix('Noise_cancelling Headphone data.xlsx','Sheet','Data','Range','Q2:Q40');

C_freq = readmatrix('Cir_output.csv','Range','A2:A2245');
C_mag_dB = readmatrix('Cir_output.csv','Range','B2:B2245');
C_phase = readmatrix('Cir_output.csv','Range','D2:D2245');

figure;
semilogx(C_freq,C_mag_dB,'k');
xlabel('Frequency (in Hz)');
ylabel('Gain (in dB)');
grid on;
%legend('Only G(s)','After adding C(s)');
title('Magnitude Bode plot of simulated circuit');

figure;
semilogx(C_freq,C_phase,'k');
xlabel('Frequency (in Hz)');
ylabel('Phase (in degree)');
grid on;
title('Phase Bode plot of simulated circuit');
%legend('Only G(s)','After adding C(s)');

Section1 = 1

%% --------------Computing gain and phase of CG and drawing its Bode plot--------------

CG_gain_dB = zeros(length(G_mag_dB));
CG_phase = zeros(length(G_phase));

for i=1:length(CG_gain_dB)
f = G_freq(i);
Cf_gain = interp1(C_freq,C_mag_dB,f);
Cf_phase = interp1(C_freq,C_phase,f);
CG_gain_dB(i) = G_mag_dB(i)+Cf_gain;
CG_phase(i) = G_phase(i)+Cf_phase;
end

% Freq = G_freq(1):10:G_freq(end);
% Freq = Freq';
% CG_gain1 = spline(G_freq,CG_gain,Freq);
% CG_phase1 = spline(G_freq,CG_phase,Freq);

figure;
semilogx(G_freq,G_mag_dB,'k--');
hold on;
semilogx(G_freq,CG_gain_dB,'k');
xlabel('Frequency (in Hz)');
ylabel('Gain (in dB)');
grid on;
legend('Only G(s)','After adding C(s)');
title('Magnitude Bode plots of G(s) and C(s)G(s)');


figure;
semilogx(G_freq,G_phase,'k--');
hold on;
semilogx(G_freq,CG_phase,'k');
xlabel('Frequency (in Hz)');
ylabel('Phase (in degree)');
grid on;
title('Phase Bode plots of G(s) and C(s)G(s)');
legend('Only G(s)','After adding C(s)');

Section2 = 2

%% ------------------G and CG-------------------------------

G = 10.^(G_mag_dB/20).*exp(1i*G_phase);
CG = 10.^(CG_gain_dB/20).*exp(1i*CG_phase);

Section3 = 3

%% ------------------B(s)/N(s)-------------------------------

B_N_tf = G./(1+CG);
B_N_before = G./(1+G);

mag = abs(B_N_tf);
phase = angle(B_N_tf);
B_N_gain_dB = 20*log10(mag(:,1));
B_N_phase = 180/pi*phase(:,1);

mag_before = abs(B_N_before);
phase_before = angle(B_N_before);
B_N_gain_dB_before = 20*log10(mag_before(:,1));
B_N_phase_before = 180/pi*phase_before(:,1);

figure;
semilogx(G_freq,B_N_gain_dB_before,'k--');
hold on;
semilogx(G_freq,B_N_gain_dB,'k');
xlabel('Frequency (in Hz)');
ylabel('Gain (in dB)');
grid on;
legend('Before adding controller','After adding controller');
title('Magnitude Bode plots of B(s)/N(s)');

Section4 = 4

%% ------------------B(s)/Vi(s)-------------------------------

B_Vi_tf = CG./(1+CG);
mag = abs(B_Vi_tf);
phase = angle(B_Vi_tf);
B_Vi_gain_dB = 20*log10(mag(:,1));
B_Vi_phase = 180/pi*phase(:,1);

figure;
semilogx(G_freq,B_Vi_gain_dB,'k');
xlabel('Frequency (in Hz)');
ylabel('Gain (in dB)');
grid on;
title('Magnitude Bode plot of B(s)/V_i(s)');

Section5 = 5