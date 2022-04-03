% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

clear

% sd = daq("ni");

model_id = ""; % id do modelo
channel = 0; % canal/canais a ler
range = 10; % range do modelo

fa = 10000; % frequência de amostragem
n_samples = 10000; % num de amostras

res_espet = fa / n_samples; % resolução espetral
res_temp = 1 / fa; % resolução temporal

t = res_temp * (0:n_samples);

% =========== Aquisição de sinais =========== %

% addinput(d, model_id, channel, "Voltage");

% data = read(d, m, "OutputFormat", "Matrix");

% data = data';

data = zeros(10, n_samples + 1);
tf = zeros(10, n_samples + 1);
tf_uni = zeros(10, n_samples / 2 + 1);

for i = 1:10
    data(i,:) = 0.8 * cos(2*pi*1000*t) + 0.010*randn(size(t));

    tf(i,:) = fft(data(i,:));
    tf(i,:) = abs(tf(i,:)) / n_samples; % módulo da transformada de fourier
    tf_uni(i,:) = tf(i,1:n_samples / 2 + 1); 
    tf_uni(i,2:end-1) = 2*tf_uni(i,2:end-1); % transformada de fourier unilateral
end

tf_uni_av = mean(tf_uni);

data_av = mean(data);

% =========== Valor Eficaz =========== %

valef = sqrt(mean(data_av.^2));

noise_ef = mean(sqrt(mean((data(1:10, :) - data_av).^2)));

% =========== Frequência =========== %

% z(1,:) = cos(2*pi*2*t);
% z(2,:) = cos(2*pi*3*t);
% z(3,:) = cos(2*pi*5*t);

% plot(t,z)

[f_fund_ef, f_fund] = max(mean(tf_uni) / sqrt(2)); % frequência fundamental e respetivo valor eficaz

f_est = (sum(tf_uni_av(f_fund-3:f_fund+3) .* (f_fund-3 - 1:f_fund+3 - 1) .* res_espet)) / ...
    sum(tf_uni_av(f_fund-3:f_fund+3)); 

%f_fund = (f_fund - 1) * res_espet;

% =========== Valor médio =========== %

% valor_medio = mean (data_av)

% =========== Espectro de potência =========== %

power_spect(1:10,:) = 20 * log10(tf_uni(1:10,:));

power_spect_den(1:10,:) = 10 * log10((tf_uni(1:10,:).^2) ./ res_espet);

% =========== THD =========== %

% harm_ef = zeros(1, (fa / 2) / (f_fund - 1) - 1);

% for i = 2:(fa / 2) / (f_fund - 1) % número máximo de harmónicas possível de amostrar
%     harm_ef(i-1) = tf_uni(1 + (f_fund - 1) * i) / sqrt(2); % valor eficaz das harmónicas
% end

% sum_harm_amp= sum(harm_ef); % soma dos valores eficazes das harmónicas

% thd_dB = 10 * log10(sum_harm_amp / f_fund_amp); % THD em dB

% =========== Diferença de fase =========== %

% u1 = cos(2*pi*50*t); % sinal 1
% u2 = cos(2*pi*50*t + pi/6); % sinal 2

% U1 = fft(u1); % transformada de fourier de sinal 1
% U2 = fft(u2); % transformada de fourier de sinal 2

% [~, index_1] = max(abs(U1)); % indíce da frequência fundamental do sinal 1
% [~, index_2] = max(abs(U2)); % indíce da frequência fundamental do sinal 2

% delta_phi = angle(U1(index_1)) - angle(U2(index_2)); % cálculo da diferença de fases

f = res_espet * (0:n_samples/2);

noise_rms = sqrt((fa / 2) * 10 .^ (mean(mean(power_spect_den)) ./ 10));

plot(f, mean(power_spect));

%figure(1);

%subplot(2, 1, 1);
%plot(t(1:100), data(1, 1:100));
%title("Dados adquiridos", "Frequência: " + f_fund + " | Valor eficaz do ruído: " + noise_ef +...
%    " | Valor eficaz: " + valef + " | Número de amostras: " + n_samples +  ...
%    " | Frequência de amostragem: " + fa + " | Alcance: " + range);

%subplot(2, 1, 2);
%plot(f, power_spect);

