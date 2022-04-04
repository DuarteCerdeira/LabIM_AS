% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

clear

model_id = ""; % id do modelo
channel = 0; % canal/canais a ler
range = 10; % range do modelo

fa = 10000; % frequência de amostragem
n_samples = 100000; % num de amostras
n_measures = 50;

res_espet = fa / n_samples; % resolução espetral
res_temp = 1 / fa; % resolução temporal

t = res_temp * (0:n_samples);

% =========== Aquisição de sinais =========== %

% sd = daq("ni");

% addinput(d, model_id, channel, "Voltage");

% data = read(d, m, "OutputFormat", "Matrix");

% data = data';

data = zeros(10, n_samples + 1);
tf = zeros(10, n_samples + 1);

for i = 1:n_measures
    data(i,:) = 0.8 * cos(2*pi*100*t) + 0.1*randn(size(t));    
    tf(i,:) = fft(data(i,:));
end

data_av = mean(data);

tf = abs(tf) / n_samples; % módulo da transformada de fourier
tf_av = mean(tf);
tf_uni = tf_av(1:n_samples / 2 + 1); 
tf_uni(2:end-1) = 2*tf_uni(2:end-1); % transformada de fourier unilateral

% =========== Frequência =========== %

[f_fund_ef, index] = max(tf_uni / sqrt(2)); % frequência fundamental e respetivo valor eficaz

f_est = (sum(tf_uni(index-3:index+3) .* ((index-1) - 3:(index-1) + 3) .* res_espet)) / ...
    sum(tf_uni(index-3:index+3)); % estimativa da frequência em caso de espalhamento espetral

f_fund = (index - 1) * res_espet;

% =========== Espectro de potência =========== %

power_spect = 20 * log10(tf_uni); % espetro de potência unilateral
power_spect_den = 10 * log10((tf_uni.^2) ./ res_espet); % densidade espetral de potência

power_spect_noise = zeros(1,n_samples/2);

for i = 1:n_samples / 2
    if i == index
        continue
    else
        power_spect_noise(i) = (tf_uni(i) ./ sqrt(2)) .^2; % espetro de potência do ruído
    end
end

power_spect_noise_den = power_spect_noise / res_espet; % densidade espetral de potência do ruído

% =========== Valor Eficaz =========== %

valef = sqrt(mean(data_av.^2));

noise_rms = sqrt((fa / 2) * mean(power_spect_noise_den));

noise_floor = 10 * log10( (noise_rms^2) / (n_samples / 2) );

% =========== PLOT =========== %

f = res_espet * (0:n_samples/2);

plot(f, power_spect);

%figure(1);

%subplot(2, 1, 1);
%plot(t(1:100), data(1, 1:100));
%title("Dados adquiridos", "Frequência: " + f_fund + " | Valor eficaz do ruído: " + noise_ef +...
%    " | Valor eficaz: " + valef + " | Número de amostras: " + n_samples +  ...
%    " | Frequência de amostragem: " + fa + " | Alcance: " + range);

%subplot(2, 1, 2);
%plot(f, power_spect);
