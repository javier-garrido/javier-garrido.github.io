close all
clear all
clc
set(0, 'DefaultLineLineWidth', 1.25, 'DefaultAxesTitleFontWeight', 'bold','DefaultAxesFontSize', 11);

%% Cargar el workspace

% [num, txt, raw] = xlsread('MEDIDAS_VNA.xlsx');
% VNArawS11_2021 = cell2table(raw(2:end, :));
% VNArawS11_2021.Properties.VariableNames = raw(1, :);
% save VNArawS11_2021.mat VNArawS11_2021
load 'VNArawS11_2021.mat'

% Temperatura de las muestras (ºC)
Temp.VNA.Air = 26.2;
Temp.VNA.DistWater = 25.9;
Temp.VNA.Acetone = 26.2;
Temp.VNA.Isopropanol = 26.2;
Temp.VNA.Methanol = 25.9;
Temp.VNA.EthyleneGlycol = 26.5;

% Gráficas del S11
col = javicolormap(12);
figure
hold on
p11 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealAir, 'Color', col(1,:), 'DisplayName', 'Real Aire', 'LineWidth', 2)
p12 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagAir, 'Color', col(2,:), 'DisplayName', 'Imag Aire', 'LineWidth', 2)
p21 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealIsopropanol, 'Color', col(3,:), 'DisplayName', 'Real Isopropanol', 'LineWidth', 2)
p22 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagIsopropanol, 'Color', col(4,:), 'DisplayName', 'Imag Isopropanol', 'LineWidth', 2)
p31 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealAcetone, 'Color', col(5,:), 'DisplayName', 'Real Acetona', 'LineWidth', 2)
p32 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagAcetone, 'Color', col(6,:), 'DisplayName', 'Imag Acetona', 'LineWidth', 2)
p41 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealMethanol, 'Color', col(7,:), 'DisplayName', 'Real Metanol', 'LineWidth', 2)
p42 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagMethanol, 'Color', col(8,:), 'DisplayName', 'Imag Metanol', 'LineWidth', 2)
p51 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealEthyleneglycol, 'Color', col(9,:), 'DisplayName', 'Real Etilenglicol', 'LineWidth', 2)
p52 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagEthyleneglycol, 'Color', col(10,:), 'DisplayName', 'Imag Etilenglicol', 'LineWidth', 2)
p61 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.RealDistWater, 'Color', col(11,:), 'DisplayName', 'Real Agua destilada', 'LineWidth', 2)
p62 = plot(VNArawS11_2021.Frequency, VNArawS11_2021.ImagDistWater, 'Color', col(12,:), 'DisplayName', 'Imag Agua destilada', 'LineWidth', 2)
legend([p11 p12 p21 p22 p31 p32 p41 p42 p51 p52 p61 p62], {'Real Aire', 'Imag Aire', 'Real Isopropanol', 'Imag Isopropanol', 'Real Acetona', 'Imag Acetona', 'Real Metanol', 'Imag Metanol', 'Real Etilenglicol', 'Imag Etilenglicol', 'Real Agua destilada', 'Imag Agua destilada'})
ylabel('S11'), xlabel('Frecuencia (Hz)'), set(gca, 'FontSize',13)

%% Interpolación para normalizar las frecuencias de los datos del S11 del VNA

Frequency = VNArawS11_2021.Frequency';
VNAS11 = table(Frequency', 'VariableNames', {'Frequency'});

% Aire
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealAir);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realAir = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagAir);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagAir = fitresult(VNAS11.Frequency);
VNAS11.Air = complex(realAir, imagAir);

% Isopropanol
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealIsopropanol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realIsopropanol = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagIsopropanol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagIsopropanol = fitresult(VNAS11.Frequency);
VNAS11.Isopropanol = complex(realIsopropanol, imagIsopropanol);

% Acetona
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealAcetone);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realAcetone = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagAcetone);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagAcetone = fitresult(VNAS11.Frequency);
VNAS11.Acetone = complex(realAcetone, imagAcetone);

% Metanol
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealMethanol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realMethanol = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagMethanol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagMethanol = fitresult(VNAS11.Frequency);
VNAS11.Methanol = complex(realMethanol, imagMethanol);

% Etilenglicol
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealEthyleneglycol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realEthyleneGlycol = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagEthyleneglycol);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagEthyleneGlycol = fitresult(VNAS11.Frequency);
VNAS11.EthyleneGlycol = complex(realEthyleneGlycol, imagEthyleneGlycol);

% Agua destilada
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.RealDistWater);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
realDistWater = fitresult(VNAS11.Frequency);
[xData, yData] = prepareCurveData(VNArawS11_2021.Frequency, VNArawS11_2021.ImagDistWater);
ft = 'pchipinterp';
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
imagDistWater = fitresult(VNAS11.Frequency);
VNAS11.DistWater = complex(realDistWater, imagDistWater);

clear xData yData fitresult ft gof realAir realIsopropanol realAcetone realDistWater realEthanol realMethanol imagAir imagIsopropanol imagAcetone imagDistWater imagEthanol imagMethanol realEthyleneGlycol imagEthyleneGlycol

%% Modelos para definir la permitividad de los materiales 
% En la bibliografía se proporcionan datos a temperaturas redondeadas.
% Necesitamos calcular los valores correspondientes para la temperatura
% de nuestras muestras, por lo que aplicaremos una interpolación lineal
% en cada caso.

% Aire
Perm_Air.VNA = ones(1, length(VNAS11.Frequency));


% Agua destilada (Debye, con datos de Kaatze (1989); interpolación entre 25ºC y 30ºC)
Eps_0_DistWater.VNA = interp1([25,30], [78.36,76.56], Temp.VNA.DistWater); 
RelaxT_DistWater.VNA = interp1([25,30], [8.27*10^-12,7.28*10^-12], Temp.VNA.DistWater);
Eps_inf_DistWater.VNA = interp1([25,30], [5.2,5.2], Temp.VNA.DistWater); 

for ii=1:length(VNAS11.Frequency)
    Perm_DistWater.VNA(ii) = Debye(Eps_0_DistWater.VNA, Eps_inf_DistWater.VNA, RelaxT_DistWater.VNA,VNAS11.Frequency(ii));
end


% Acetona (Cole-Cole, con datos de Buckley and Maryott (1958) (página 8); interpolación entre 20ºC y 40ºC)
Eps_0_Acetone.VNA = interp1([1,20,40], [23.29,21.2,19.29], Temp.VNA.Acetone);
Eps_inf_Acetone.VNA = interp1([1,20,40], [1.93,1.9,0.87], Temp.VNA.Acetone); 
RelaxT_Acetone.VNA = interp1([1,20,40], [0.75,0.63,0.52], Temp.VNA.Acetone)*10^-2/(2*pi*physconst('LightSpeed')); 

for ii=1:length(VNAS11.Frequency)
    Perm_Acetone.VNA(ii) = ColeCole(Eps_0_Acetone.VNA, Eps_inf_Acetone.VNA, RelaxT_Acetone.VNA, 1, VNAS11.Frequency(ii));
end


% Isopropanol (Double Debye, con datos de Gregory and Clarke (2012); interpolación entre 25ºC y 30ºC)
Eps_s_Isopropanol.VNA = interp1([25,30], [19.3,18.5], Temp.VNA.Isopropanol);
fr1_Isopropanol.VNA = interp1([25,30], [0.443E9,0.557E9], Temp.VNA.Isopropanol);
Eps_h_Isopropanol.VNA = interp1([25,30], [3.551,3.538], Temp.VNA.Isopropanol);
fr2_Isopropanol.VNA = interp1([25,30], [5.999E9,6.655E9], Temp.VNA.Isopropanol);
Eps_inf_Isopropanol.VNA = interp1([25,30], [3.065,3.056], Temp.VNA.Isopropanol);

for ii=1:length(VNAS11.Frequency)
    Perm_Isopropanol.VNA(ii) = Eps_inf_Isopropanol.VNA + (Eps_s_Isopropanol.VNA-Eps_h_Isopropanol.VNA) / (1+j*VNAS11.Frequency(ii)/fr1_Isopropanol.VNA) + (Eps_h_Isopropanol.VNA-Eps_inf_Isopropanol.VNA) / (1+j*VNAS11.Frequency(ii)/fr2_Isopropanol.VNA);
end


% Metanol (Single Debye, con datos de Gregory and Clarke (2012); interpolación entre 25ºC y 30ºC)
Eps_s_Methanol.VNA = interp1([25,30], [32.66,31.69], Temp.VNA.Methanol);
fr_Methanol.VNA = interp1([25,30], [3.141E9,3.49E9], Temp.VNA.Methanol);
Eps_inf_Methanol.VNA = interp1([25,30], [5.563,5.45], Temp.VNA.Methanol);

for ii=1:length(VNAS11.Frequency)
    Perm_Methanol.VNA(ii) = Eps_inf_Methanol.VNA + (Eps_s_Methanol.VNA-Eps_inf_Methanol.VNA) / (1+j*VNAS11.Frequency(ii) / fr_Methanol.VNA);
end

% Etilenglicol (Davidson-Cole, con datos de Gregory and Clarke (2012); interpolación entre 25ºC y 30ºC) 
Eps_0_EthyleneGlycol.VNA = interp1([25,30], [40.75,39.66], Temp.VNA.EthyleneGlycol); 
Eps_inf_EthyleneGlycol.VNA = interp1([25,30], [4.7,4.712], Temp.VNA.EthyleneGlycol); 
Relaxfreq_EthyleneGlycol.VNA = interp1([25,30], [1.19E9,1.451E9], Temp.VNA.EthyleneGlycol);
Beta.VNA = interp1([25,30], [0.859,0.864], Temp.VNA.EthyleneGlycol); 

for ii=1:length(VNAS11.Frequency)
    Perm_EthyleneGlycol.VNA(ii) = ColeDavidson(Eps_0_EthyleneGlycol.VNA, Eps_inf_EthyleneGlycol.VNA, Relaxfreq_EthyleneGlycol.VNA,Beta.VNA, VNAS11.Frequency(ii));
end

% Comparativa de las permitividades teóricas para la temperatura
% de las muestras
col = javicolormap(7);

figure
hold on
p1 = plot(VNAS11.Frequency, real(Perm_Air.VNA), 'Color',col(1,:), 'DisplayName', 'Aire', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_Air.VNA), 10, '*', 'MarkerEdgeColor', col(1,:))
p2 = plot(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 'Color', col(2,:), 'DisplayName', 'Isopropanol', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 10, '*', 'MarkerEdgeColor', col(2,:))
p3 = plot(VNAS11.Frequency, real(Perm_Acetone.VNA), 'Color', col(3,:), 'DisplayName', 'Acetona', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_Acetone.VNA), 10, '*', 'MarkerEdgeColor', col(3,:))
p4 = plot(VNAS11.Frequency, real(Perm_Methanol.VNA), 'Color', col(5,:), 'DisplayName', 'Metanol', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_Methanol.VNA), 10, '*', 'MarkerEdgeColor', col(5,:))
p5 = plot(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 'Color', col(6,:), 'DisplayName', 'Etilenglicol', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 10, '*', 'MarkerEdgeColor', col(6,:))
p6 = plot(VNAS11.Frequency, real(Perm_DistWater.VNA), 'Color', col(7,:), 'DisplayName', 'Agua destilada', 'LineWidth', 2)
scatter(VNAS11.Frequency, real(Perm_DistWater.VNA), 10, '*', 'MarkerEdgeColor', col(7,:))
legend([p1 p2 p3 p4 p5 p6], {'Aire', 'Isopropanol', 'Acetona', 'Metanol', 'Etilenglicol', 'Agua destilada'})
ylabel('\epsilon'), xlabel('Frecuencia (Hz)'), set(gca, 'FontSize',13)

figure
hold on
p1 = plot(VNAS11.Frequency, -imag(Perm_Air.VNA), 'Color', col(1,:), 'DisplayName', 'Aire', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_Air.VNA), 10, '*', 'MarkerEdgeColor', col(1,:))
p2 = plot(VNAS11.Frequency, -imag(Perm_Isopropanol.VNA), 'Color', col(2,:), 'DisplayName', 'Isopropanol', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_Isopropanol.VNA), 10, '*', 'MarkerEdgeColor', col(2,:))
p3 = plot(VNAS11.Frequency, -imag(Perm_Acetone.VNA), 'Color', col(3,:), 'DisplayName', 'Acetona', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_Acetone.VNA), 10, '*', 'MarkerEdgeColor', col(3,:))
p4 = plot(VNAS11.Frequency, -imag(Perm_Methanol.VNA), 'Color', col(5,:), 'DisplayName', 'Metanol', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_Methanol.VNA), 10, '*', 'MarkerEdgeColor', col(5,:))
p5 = plot(VNAS11.Frequency, -imag(Perm_EthyleneGlycol.VNA), 'Color', col(6,:), 'DisplayName', 'Etilenglicol', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_EthyleneGlycol.VNA), 10, '*', 'MarkerEdgeColor', col(6,:))
p6 = plot(VNAS11.Frequency, -imag(Perm_DistWater.VNA), 'Color', col(7,:), 'DisplayName', 'Agua destilada', 'LineWidth', 2)
scatter(VNAS11.Frequency, -imag(Perm_DistWater.VNA), 10, '*', 'MarkerEdgeColor', col(7,:))
legend([p1 p2 p3 p4 p5 p6], {'Aire', 'Isopropanol', 'Acetona', 'Metanol', 'Etilenglicol', 'Agua destilada'})
ylabel('\epsilon'), xlabel('Frecuencia (Hz)'), set(gca, 'FontSize',13)

%% Permitividad a partir del S11 del VNA
% Seguimos el procedimiento descrito en "Wagner et al. Numerical 3-D FEM
% and Experimental Analysis of the Open-Ended Coaxial Line Technique for
% Microwave Dielectric Spectroscopy on Soil", 2014.

% -> Calibración Open-Water-Liquid (OWL)


% Calibración con Open(Aire), Water(Agua destilada) y Liquid (Metanol).
for ii=1:length(VNAS11.Frequency)
    M(:,:,ii) = [VNAS11.Air(ii) -1 -Perm_Air.VNA(ii);VNAS11.DistWater(ii) -1 -Perm_DistWater.VNA(ii);VNAS11.Methanol(ii) -1 -Perm_Methanol.VNA(ii)];
    e(:,ii) = [-Perm_Air.VNA(ii)*VNAS11.Air(ii); -Perm_DistWater.VNA(ii)*VNAS11.DistWater(ii); -Perm_Methanol.VNA(ii)*VNAS11.Methanol(ii)];
    c(:,ii) = inv(M(:,:,ii))*e(:,ii);
end

for ii=1:length(VNAS11.Frequency)
    Estim_VNA_Perm_Air(ii) = transpose((c(1,ii)*VNAS11.Air(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Air(ii)));
    Estim_VNA_Perm_Isopropanol(ii) = transpose((c(1,ii)*VNAS11.Isopropanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Isopropanol(ii)));
    Estim_VNA_Perm_Acetone(ii) = transpose((c(1,ii)*VNAS11.Acetone(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Acetone(ii)));
    Estim_VNA_Perm_Methanol(ii) = transpose((c(1,ii)*VNAS11.Methanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Methanol(ii)));
    Estim_VNA_Perm_EthyleneGlycol(ii) = transpose((c(1,ii)*VNAS11.EthyleneGlycol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.EthyleneGlycol(ii)));
    Estim_VNA_Perm_DistWater(ii) = transpose((c(1,ii)*VNAS11.DistWater(ii)-c(2,ii)) / (c(3,ii)-VNAS11.DistWater(ii)));    
end

SE_Methanol_VNA = sum([abs(real(Estim_VNA_Perm_Air) - real(Perm_Air.VNA)),abs(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)),...
    abs(real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)),...
    abs(real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)), abs(real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)),...
    abs(real(Estim_VNA_Perm_DistWater) - real(Perm_DistWater.VNA))]); 

RMSE_Methanol_VNA = sqrt((sum([(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)).^2,...
    (real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)).^2,...
    (real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)).^2]))/(3*length(Estim_VNA_Perm_Air))); 

figure
hold on
plot(VNAS11.Frequency,real(Perm_Air.VNA), 'Color', col(1,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Air), 10, '*', 'MarkerEdgeColor', col(1,:))
plot(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 'Color', col(2,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Isopropanol), 10, '*', 'MarkerEdgeColor', col(2,:))
plot(VNAS11.Frequency, real(Perm_Acetone.VNA), 'Color', col(3,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Acetone), 10, '*', 'MarkerEdgeColor', col(3,:))
plot(VNAS11.Frequency, real(Perm_Methanol.VNA), 'Color', col(5,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Methanol), 10, '*', 'MarkerEdgeColor', col(5,:))
plot(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 'Color', col(6,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_EthyleneGlycol), 10, '*', 'MarkerEdgeColor', col(6,:))
plot(VNAS11.Frequency, real(Perm_DistWater.VNA), 'Color', col(7,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_DistWater), 10, '*', 'MarkerEdgeColor', col(7,:))
ylim([0 80])
%xlim([0 900e6])
title('Calibración OWL con Metanol')
ylabel('Re(\epsilon)')
xlabel('Frecuencia (Hz)')
legend('Aire modelo', 'Aire estimada', 'Isopropanol modelo', 'Isopropanol estimada', 'Acetona modelo', 'Acetona estimada', 'Metanol modelo', 'Metanol estimada', 'Etilenglicol modelo', 'Etilenglicol estimada', 'Agua destilada modelo', 'Agua destilada estimada')


% Calibración con Open(Aire), Water(Agua destilada) y Liquid (Isopropanol)
for ii=1:length(VNAS11.Frequency)
    M(:,:,ii) = [VNAS11.Air(ii) -1 -Perm_Air.VNA(ii);VNAS11.DistWater(ii) -1 -Perm_DistWater.VNA(ii);VNAS11.Isopropanol(ii) -1 -Perm_Isopropanol.VNA(ii)];
    e(:,ii) = [-Perm_Air.VNA(ii)*VNAS11.Air(ii); -Perm_DistWater.VNA(ii)*VNAS11.DistWater(ii); -Perm_Isopropanol.VNA(ii)*VNAS11.Isopropanol(ii);];
    c(:,ii) = inv(M(:,:,ii))*e(:,ii);
end

for ii=1:length(VNAS11.Frequency)
    Estim_VNA_Perm_Air(ii) = transpose((c(1,ii)*VNAS11.Air(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Air(ii)));
    Estim_VNA_Perm_Isopropanol(ii) = transpose((c(1,ii)*VNAS11.Isopropanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Isopropanol(ii)));
    Estim_VNA_Perm_Acetone(ii) = transpose((c(1,ii)*VNAS11.Acetone(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Acetone(ii)));
    Estim_VNA_Perm_Methanol(ii) = transpose((c(1,ii)*VNAS11.Methanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Methanol(ii)));
    Estim_VNA_Perm_EthyleneGlycol(ii) = transpose((c(1,ii)*VNAS11.EthyleneGlycol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.EthyleneGlycol(ii)));
    Estim_VNA_Perm_DistWater(ii) = transpose((c(1,ii)*VNAS11.DistWater(ii)-c(2,ii)) / (c(3,ii)-VNAS11.DistWater(ii)));    
end
SE_Isopropanol_VNA = sum([abs(real(Estim_VNA_Perm_Air) - real(Perm_Air.VNA)),abs(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)),...
    abs(real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)),...
    abs(real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)), abs(real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)),...
    abs(real(Estim_VNA_Perm_DistWater) - real(Perm_DistWater.VNA))]); 
RMSE_Isopropanol_VNA = sqrt((sum([(real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)).^2,...
    (real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)).^2,...
    (real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)).^2]))/(3*length(Estim_VNA_Perm_Air))); 

figure
hold on
plot(VNAS11.Frequency, real(Perm_Air.VNA), 'Color', col(1,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Air), 10, '*', 'MarkerEdgeColor', col(1,:))
plot(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 'Color', col(2,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Isopropanol), 10, '*', 'MarkerEdgeColor', col(2,:))
plot(VNAS11.Frequency, real(Perm_Acetone.VNA), 'Color', col(3,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Acetone), 10, '*', 'MarkerEdgeColor', col(3,:))
plot(VNAS11.Frequency, real(Perm_Methanol.VNA), 'Color', col(5,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Methanol), 10, '*', 'MarkerEdgeColor', col(5,:))
plot(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 'Color', col(6,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_EthyleneGlycol), 10, '*', 'MarkerEdgeColor', col(6,:))
plot(VNAS11.Frequency, real(Perm_DistWater.VNA), 'Color', col(7,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_DistWater), 10, '*', 'MarkerEdgeColor', col(7,:))
ylim([0 80])
%xlim([0 900e6])
title('Calibración OWL con Isopropanol')
ylabel('Re(\epsilon)')
xlabel('Frecuencia (Hz)')
legend('Aire modelo', 'Aire estimada', 'Isopropanol modelo', 'Isopropanol estimada', 'Acetona modelo', 'Acetona estimada', 'Metanol modelo', 'Metanol estimada', 'Etilenglicol modelo', 'Etilenglicol estimada', 'Agua destilada modelo', 'Agua destilada estimada')


% Calibración con Open(Aire), Water(Agua destilada) y Liquid (Acetona)
% Este es el mejor (error más bajo)
for ii=1:length(VNAS11.Frequency)
    M(:,:,ii) = [VNAS11.Air(ii) -1 -Perm_Air.VNA(ii);VNAS11.DistWater(ii) -1 -Perm_DistWater.VNA(ii);VNAS11.Acetone(ii) -1 -Perm_Acetone.VNA(ii)];
    e(:,ii) = [-Perm_Air.VNA(ii)*VNAS11.Air(ii); -Perm_DistWater.VNA(ii)*VNAS11.DistWater(ii); -Perm_Acetone.VNA(ii)*VNAS11.Acetone(ii);];
    c(:,ii) = inv(M(:,:,ii))*e(:,ii);
end

for ii=1:length(VNAS11.Frequency)
    Estim_VNA_Perm_Air(ii) = transpose((c(1,ii)*VNAS11.Air(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Air(ii)));
    Estim_VNA_Perm_Isopropanol(ii) = transpose((c(1,ii)*VNAS11.Isopropanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Isopropanol(ii)));
    Estim_VNA_Perm_Acetone(ii) = transpose((c(1,ii)*VNAS11.Acetone(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Acetone(ii)));
    Estim_VNA_Perm_Methanol(ii) = transpose((c(1,ii)*VNAS11.Methanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Methanol(ii)));
    Estim_VNA_Perm_EthyleneGlycol(ii) = transpose((c(1,ii)*VNAS11.EthyleneGlycol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.EthyleneGlycol(ii)));
    Estim_VNA_Perm_DistWater(ii) = transpose((c(1,ii)*VNAS11.DistWater(ii)-c(2,ii)) / (c(3,ii)-VNAS11.DistWater(ii)));    
end
SE_Acetone_VNA = sum([abs(real(Estim_VNA_Perm_Air) - real(Perm_Air.VNA)),abs(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)),...
    abs(real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)),...
    abs(real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)),abs(real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)),...
    abs(real(Estim_VNA_Perm_DistWater) - real(Perm_DistWater.VNA))]); 
RMSE_Acetone_VNA = sqrt((sum([(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)).^2,...
    (real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)).^2,...
    (real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)).^2]))/(3*length(Estim_VNA_Perm_Air))); 

figure
hold on
plot(VNAS11.Frequency, real(Perm_Air.VNA), 'Color', col(1,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Air), 10, '*', 'MarkerEdgeColor', col(1,:))
plot(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 'Color', col(2,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Isopropanol), 10, '*', 'MarkerEdgeColor', col(2,:))
plot(VNAS11.Frequency, real(Perm_Acetone.VNA), 'Color', col(3,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Acetone), 10, '*', 'MarkerEdgeColor', col(3,:))
plot(VNAS11.Frequency, real(Perm_Methanol.VNA), 'Color', col(5,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Methanol), 10, '*', 'MarkerEdgeColor', col(5,:))
plot(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 'Color', col(6,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_EthyleneGlycol), 10, '*', 'MarkerEdgeColor', col(6,:))
plot(VNAS11.Frequency, real(Perm_DistWater.VNA), 'Color', col(7,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_DistWater), 10, '*', 'MarkerEdgeColor', col(7,:))
ylim([0 80])
%xlim([0 900e6])
title('Calibración OWL con Acetona')
ylabel('Re(\epsilon)')
xlabel('Frecuencia (Hz)')
legend('Aire modelo', 'Aire estimada', 'Isopropanol modelo', 'Isopropanol estimada', 'Acetona modelo', 'Acetona estimada', 'Metanol modelo', 'Metanol estimada', 'Etilenglicol modelo', 'Etilenglicol estimada', 'Agua destilada modelo', 'Agua destilada estimada')


% Calibración con Open(Aire), Water(Agua destilada) y Liquid (Etilenglicol)
for ii=1:length(VNAS11.Frequency)
    M(:,:,ii) = [VNAS11.Air(ii) -1 -Perm_Air.VNA(ii);VNAS11.DistWater(ii) -1 -Perm_DistWater.VNA(ii);VNAS11.EthyleneGlycol(ii) -1 -Perm_EthyleneGlycol.VNA(ii)];
    e(:,ii) = [-Perm_Air.VNA(ii)*VNAS11.Air(ii); -Perm_DistWater.VNA(ii)*VNAS11.DistWater(ii); -Perm_EthyleneGlycol.VNA(ii)*VNAS11.EthyleneGlycol(ii);];
    c(:,ii) = inv(M(:,:,ii))*e(:,ii);
end

for ii=1:length(VNAS11.Frequency)
    Estim_VNA_Perm_Air(ii) = transpose((c(1,ii)*VNAS11.Air(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Air(ii)));
    Estim_VNA_Perm_Isopropanol(ii) = transpose((c(1,ii)*VNAS11.Isopropanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Isopropanol(ii)));
    Estim_VNA_Perm_Acetone(ii) = transpose((c(1,ii)*VNAS11.Acetone(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Acetone(ii)));
    Estim_VNA_Perm_Methanol(ii) = transpose((c(1,ii)*VNAS11.Methanol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.Methanol(ii)));
    Estim_VNA_Perm_EthyleneGlycol(ii) = transpose((c(1,ii)*VNAS11.EthyleneGlycol(ii)-c(2,ii)) / (c(3,ii)-VNAS11.EthyleneGlycol(ii)));
    Estim_VNA_Perm_DistWater(ii) = transpose((c(1,ii)*VNAS11.DistWater(ii)-c(2,ii)) / (c(3,ii)-VNAS11.DistWater(ii)));    
end
SE_EthyleneGlycol_VNA = sum([abs(real(Estim_VNA_Perm_Air) - real(Perm_Air.VNA)),abs(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)),...
    abs(real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)),...
    abs(real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)), abs(real(Estim_VNA_Perm_EthyleneGlycol) - real(Perm_EthyleneGlycol.VNA)),...
    abs(real(Estim_VNA_Perm_DistWater) - real(Perm_DistWater.VNA))]); 
RMSE_EthyleneGlycol_VNA = sqrt((sum([(real(Estim_VNA_Perm_Isopropanol) - real(Perm_Isopropanol.VNA)).^2,...
    (real(Estim_VNA_Perm_Acetone) - real(Perm_Acetone.VNA)).^2,...
    (real(Estim_VNA_Perm_Methanol) - real(Perm_Methanol.VNA)).^2]))/(3*length(Estim_VNA_Perm_Air))); 

figure
hold on
plot(VNAS11.Frequency, real(Perm_Air.VNA), 'Color', col(1,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Air), 10, '*', 'MarkerEdgeColor', col(1,:))
plot(VNAS11.Frequency, real(Perm_Isopropanol.VNA), 'Color', col(2,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Isopropanol), 10, '*', 'MarkerEdgeColor', col(2,:))
plot(VNAS11.Frequency, real(Perm_Acetone.VNA), 'Color', col(3,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Acetone), 10, '*', 'MarkerEdgeColor', col(3,:))
plot(VNAS11.Frequency, real(Perm_Methanol.VNA), 'Color', col(5,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_Methanol), 10, '*', 'MarkerEdgeColor', col(5,:))
plot(VNAS11.Frequency, real(Perm_EthyleneGlycol.VNA), 'Color', col(6,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_EthyleneGlycol), 10, '*', 'MarkerEdgeColor', col(6,:))
plot(VNAS11.Frequency, real(Perm_DistWater.VNA), 'Color', col(7,:))
scatter(VNAS11.Frequency, real(Estim_VNA_Perm_DistWater), 10, '*', 'MarkerEdgeColor', col(7,:))
ylim([0 80])
%xlim([0 900e6])
title('Calibración OWL con Etilenglicol')
ylabel('Re(\epsilon)')
xlabel('Frecuencia (Hz)')
legend('Aire modelo', 'Aire estimada', 'Isopropanol modelo', 'Isopropanol estimada', 'Acetona modelo', 'Acetona estimada', 'Metanol modelo', 'Metanol estimada', 'Etilenglicol modelo', 'Etilenglicol estimada', 'Agua destilada modelo', 'Agua destilada estimada')



%% Funciones
function [Debye_Perm] = Debye(eps_0, eps_inf, relaxT, freq)
% Función que calcula el modelo de Debye de medios dieléctricos
% La función ofrece como salida el valor complejo de la permitividad a partir de las entradas: 1) 'eps_0' = Permitividad estática; 2) 'eps_inf' = Permitividad en el infinito; 3) 'relaxT' = Tiempo de relajación (s); 4) 'freq' = Frecuencia (Hz)
    
    Debye_Perm = eps_inf+(eps_0 - eps_inf) / (1 + j*2*pi*freq*relaxT);
end

function [ColeCole_Perm] = ColeCole(eps_0, eps_inf, relaxT, beta, freq)
% Función que calcula el modelo de Cole-Cole de medios dieléctricos
% La función ofrece como salida el valor complejo de la permitividad a partir de las entradas: 1) 'eps_0' = Permitividad estática; 2) 'eps_inf' = Permitividad en el infinito; 3) 'relaxT' = Tiempo de relajación (s); 4) 'beta' = Parámetro que caracteriza los tiempos de la distribución de relajación; 5) 'freq' = Frecuencia (Hz)

    ColeCole_Perm = eps_inf+(eps_0 - eps_inf) / (1 + (j*2*pi*freq*relaxT)^beta);
end

function [ColeDavidson_Perm] = ColeDavidson(eps_0, eps_inf, relaxfreq, beta, freq)
% Función que calcula el modelo de Davidson-Cole de medios dieléctricos
% La función ofrece como salida el valor complejo de la permitividad a partir de las entradas: 1) 'eps_0' = Permitividad estática; 2) 'eps_inf' = Permitividad en el infinito; 3) 'relaxT' = Tiempo de relajación (s); 4) 'beta' = Parámetro que caracteriza los tiempos de la distribución de relajación; 5) 'freq' = Frecuencia (Hz)

    ColeDavidson_Perm = eps_inf+(eps_0 - eps_inf) / (1 + j*freq/relaxfreq)^beta;
end







function J = javicolormap(m)
    % Mapa de colores
    
    c = [1, 60, 244;
        239, 1, 51;
        48, 239, 36;
        194, 85, 230;
        252, 166, 26;
        29, 214, 255;
        254, 139, 255;
        41, 238, 162;
        141, 24, 249;
        247, 193, 34;
        20, 128, 253;
        253, 44, 253;
        161, 252, 55;
        96, 4, 199;
        168, 104, 0];

    if m <= length(c)
        J = zeros(m,3);
        for ii=1:m
            J(ii,:) = c(ii,:)./255;
        end
    else
        J = zeros(m,3);
        positions = linspace(0,1,length(c));
        J = customcolormap(positions,c,m);
        J = J./255;
    end
end