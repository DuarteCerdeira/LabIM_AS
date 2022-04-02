% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

d = daq("ni");

model_id = "";
channel = 0;
range = 10;

fa = 1000; % frequência de amostragem
n_samples = 0; % num de amostras

res_espet = fa / n_samples;
res_temp = 1 / fa;


addinput(d, model_id, channel, "Voltage");

data = read(d, m, "OutputFormat", "Matrix");

