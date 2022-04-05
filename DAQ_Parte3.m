% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

clear

% sd = daq("ni");

model_id = ""; % id do modelo
channel = 0; % canal/canais a ler
range = 10; % range do modelo

fa = 10000; % frequência de amostragem
n_samples = 5000; % num de amostras

res_espet = fa / n_samples; % resolução espetral
res_temp = 1 / fa; % resolução temporal

t = res_temp * (0:n_samples);

% =========== Aquisição de sinais =========== %

% addinput(d, model_id, channel, "Voltage");

% data = read(d, m, "OutputFormat", "Matrix");

% data = data';

data(1,:) = 3 * cos(2*pi*100*t + pi/4);
data(2,:) = 5 * cos(2*pi*100*t);

% =========== Valor Eficaz =========== %

valef = (sqrt(mean(data'.^2)))';

% =========== Frequência =========== %

tf(1,:) = fft(data(1,:));
tf(2,:) = fft(data(2,:));
tf_mod = abs(tf) / n_samples; % módulo da transformada de fourier
tf_uni = tf_mod(:,1:n_samples / 2 + 1); 
tf_uni(2:end-1) = 2*tf_uni(2:end-1); % transformada de fourier unilateral

[f_fund_ef, index] = max(tf_uni, [], 2); % frequência fundamental e respetivo valor eficaz

f_est(1) = (sum(tf_uni(1, index(1)-3:index(1)+3) .* ((index(1)-1) - 3:(index(1)-1) + 3) .* res_espet)) / ...
    sum(tf_uni(1, index(1)-3:index(1)+3)); % estimativa da frequência em caso de espalhamento espetral

f_est(2) = (sum(tf_uni(2, index(2)-3:index(2)+3) .* ((index(2)-1) - 3:(index(2)-1) + 3) .* res_espet)) / ...
    sum(tf_uni(2, index(2)-3:index(2)+3)); % estimativa da frequência em caso de espalhamento espetral

f_fund = (index - 1) * res_espet;

% =========== Diferença de fase =========== %

delta_phi = angle(tf(1,index(1))) - angle(tf(2,index(2))); % cálculo da diferença de fases

f = res_espet * (0:n_samples/2);

figure(1);

subplot(2, 1, 1);
plot(t(1:500), data(1, 1:500));
title("Dados adquiridos", "Frequência: " + f_fund(1) + " | Diferença de fase: " + rad2deg(delta_phi) + ...
    " | Valor eficaz: " + valef(1) + " | Número de amostras: " + n_samples +  ...
    " | Frequência de amostragem: " + fa + " | Alcance: " + range);

subplot(2, 1, 2);
plot(t(1:500), data(2, 1:500));
title("Dados adquiridos", "Frequência: " + f_fund(2) + " | Diferença de fase: " + rad2deg(delta_phi) + ...
    " | Valor eficaz: " + valef(2) + " | Número de amostras: " + n_samples +  ...
    " | Frequência de amostragem: " + fa + " | Alcance: " + range);

