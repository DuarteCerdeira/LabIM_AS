% Instrumentação e Medidas - Laboratório 5 - Aquisição de Sinais
% 96195 - Duarte Cerdeira
% Outros caralhos

d = daq("ni");

model_id = "";
channel = 0;

fa = 1000; % frequência de amostragem
m = 0; % num de amostras

addinput(d, model_id, channel, "Voltage");

data = read(d, m, "OutputFormat", "Matrix");

