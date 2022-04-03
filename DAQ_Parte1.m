% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

clear

d = daq("ni");

model_id = ""; % id do modelo
channel = 0; % canal/canais a ler
range = 10; % range do modelo

fa = 1000; % frequência de amostragem
n_samples = 5000; % num de amostras

res_espet = fa / n_samples; % resolução espetral
res_temp = 1 / fa; % resolução temporal

tf = []; % transformada de 
tf_uni = [];
power_spect = [];

addinput(d, model_id, channel, "Voltage");

data = read(d, m, "OutputFormat", "Matrix");

% =========== Valor Eficaz =========== %

media = mean((data.^2))
valef = sqrt(media)

% =========== Frequência =========== %

res_espet = 1

z(1,:) = cos(2*pi*2*t)
z(2,:) = cos(2*pi*3*t)
z(3,:) = cos(2*pi*5*t)

plot(t,z)

freqs = zeros(1, height(z))

for index = 1:height(z)
    tf = fft(z(index,:));
    tf = abs(tf / n_samples);
    tf_uni = tf(1:n_samples / 2 + 1);
    tf_uni(2:end-1) = 2*tf_uni(2:end-1);
    f = res_espet * (0:n_samples/2);
    [~,I] = max(abs(tf));
    freqs(index) = freq(I)
end

plot(f, tf_uni);
xlim([0 10])

% =========== Valor médio =========== %

valor_medio = mean (data);

% =========== Espectro =========== %

signal_espetro = tf_uni;
signal_espetro_P = signal_espetro^2

% =========== THD =========== %

[f_fund_ef, f_fund] = max(tf_uni ./ sqrt(2)); % frequência fundamental e respetivo valor eficaz

harm_ef = zeros(1, (fa / 2) / (f_fund - 1) - 1);

for i = 2:(fa / 2) / (f_fund - 1) % número máximo de harmónicas possível de amostrar
    harm_ef(i-1) = tf_uni(1 + (f_fund - 1) * i) / sqrt(2); % valor eficaz das harmónicas
end

sum_harm_amp= sum(harm_ef); % soma dos valores eficazes das harmónicas

thd_dB = 10 * log10(sum_harm_amp / f_fund_amp); % THD em dB

% =========== Diferença de fase =========== %

u1 = cos(2*pi*50*t); % sinal 1
u2 = cos(2*pi*50*t + pi/6); % sinal 2

U1 = fft(u1); % transformada de fourier de sinal 1
U2 = fft(u2); % transformada de fourier de sinal 2

[~, index_1] = max(abs(U1)); % indíce da frequência fundamental do sinal 1
[~, index_2] = max(abs(U2)); % indíce da frequência fundamental do sinal 2

delta_phi = angle(U1(index_1)) - angle(U2(index_2)); % cálculo da diferença de fases


