%% NBR 15751:2013 - Anexo B - Calculo das tensoes de passo e de toque
%% CASO 2 - novo solo estratificado (rho2 > rho1)
%% ========================================================================
%% Mesma malha e mesma metodologia do NBR_Potenciais_0A, mas com o
%% novo modelo de solo: rho_1=197,972 ohm.m (h1=2,636 m), rho_2=787,754
%% ohm.m (infinita), rho_a=524,125 ohm.m.
%% ========================================================================
%% -----------------------------------------------------------------------
%% 1) Dados reaproveitados (inalterados em relacao ao Caso 1)
%% -----------------------------------------------------------------------
If = 12.5;              % Corrente de falta fase-terra, kA
t  = 0.3;                % Tempo de atuacao do sistema de protecao, s

Df = 1.085;              % Fator de decremento
Sf = 0.8;                % Fator de distribuicao
Cp = 1.1;                % Fator de projecao (crescimento do sistema)

Im_kA = Df*Sf*If*Cp;     % Corrente de malha simetrica eficaz [19], kA
Imcd  = Im_kA*1000;      % Convertida para A (curta duracao) - usada no Anexo B
fprintf('Corrente de malha (Im) = %.4f kA = %.1f A\n', Im_kA, Imcd);

Rch  = 1000;              % Resistencia do corpo humano, ohm
Ichcd = 0.116/sqrt(t);    % Corrente de choque de curta duracao, A [6]

%% -----------------------------------------------------------------------
%% 2) Novo modelo de solo (CASO 2)
%% -----------------------------------------------------------------------
rho_s = 3000;             % Resistividade da camada superficial (brita), ohm.m
rho_1 = 197.972;           % Resistividade da 1a camada de solo, ohm.m
h1    = 2.636;              % Profundidade da 1a camada, m
rho_2 = 787.754;            % Resistividade da 2a camada de solo (infinita), ohm.m
rho_a = 524.125;            % Resistividade aparente - DADA diretamente pelo autor

a2 = 1-(rho_1/rho_s);
a3 = 2*0.1+0.106;          % a = 0,106 m (constante da NBR, eq [14])
C  = 1 - 0.106*(a2/a3);    % Fator de reducao C
fprintf('Fator de reducao C = %.4f\n', C);

Valor_permissivel_passo = (Rch + 6*rho_s*C)*Ichcd;   %#ok<*NOPTS>  % [10]
Valor_permissivel_toque = (Rch + 1.5*rho_s*C)*Ichcd;              % [11]
fprintf('Tensao de passo permissivel = %.1f V\n', Valor_permissivel_passo);
fprintf('Tensao de toque permissivel = %.1f V\n', Valor_permissivel_toque);

%% -----------------------------------------------------------------------
%% 3) Resistividade do solo a usar no calculo dos potenciais (Anexo B)
%% -----------------------------------------------------------------------
% O Anexo B (texto introdutorio) recomenda:
%   a) usar rho_1 quando rho_2 < rho_1;
%   b) usar rho_a (aparente) quando rho_2 > rho_1.
% AGORA rho_2 (787,754) > rho_1 (197,972) -> usa-se rho_a (dado direto).
if rho_2 < rho_1
    rho_calc = rho_1;
    fprintf('rho_2 < rho_1 -> usando rho_1 = %.3f ohm.m no calculo de Vt/Vp\n', rho_calc);
else
    rho_calc = rho_a;
    fprintf('rho_2 > rho_1 -> usando rho_a = %.3f ohm.m no calculo de Vt/Vp\n', rho_calc);
end

%% -----------------------------------------------------------------------
%% 4) Geometria real da malha (extraida do DXF - inalterada)
%% -----------------------------------------------------------------------
Lt_horizontais = 11241.416;   % Comprimento dos condutores horizontais, m (do DXF)
n_hastes = 48;
comp_haste = 3.0;             % m

Lt = Lt_horizontais + n_hastes*comp_haste;   % B.1.1: inclui hastes
fprintf('Lt (condutores + hastes) = %.3f m\n', Lt);

H = 0.6;                      % Profundidade da malha, m - ASSUMIDA (ver aviso no topo)
Secao_utilizada_mm2 = 35;     % Secao nominal - igual ao Caso 1 (mesmo If, t)
d = 7.5e-3;                    % Diametro real do condutor (XGSLab), m
fprintf('Diametro do condutor = %.4f mm\n', d*1000);

%% -----------------------------------------------------------------------
%% 5) Fator kii (inalterado - mesma geometria de hastes)
%% -----------------------------------------------------------------------
hastes_com_periferia = true;

%% -----------------------------------------------------------------------
%% 6) Calculo de Km, Ki, Ks e Vt/Vp nas duas direcoes da malha
%% -----------------------------------------------------------------------
D = [6.25, 12.5];
n = [25, 25];

Vt = zeros(1,2);
Vp = zeros(1,2);
Km = zeros(1,2);
Ks = zeros(1,2);
Ki = zeros(1,2);
kii = zeros(1,2);

for k = 1:2
    Dk = D(k);
    nk = n(k);

    if hastes_com_periferia
        kii(k) = 1;
    else
        kii(k) = (2*nk)^(-2/nk);         % [26]
    end

    Ki(k) = 0.656 + 0.172*nk;             % [27]

    Km(k) = (1/(2*pi)) * ( ...
                log( (Dk^2)/(16*H*d) + ((Dk+2*H)^2)/(8*Dk*d) - H/(4*d) ) ...
                + (kii(k)/sqrt(1+H)) * log(8/(pi*(2*nk-1))) );   % [25]

    Ks(k) = (1/pi) * ( 1/(2*H) + 1/(Dk+H) + (1/Dk)*(1-0.5^(nk-2)) );  % [30]

    Vt(k) = (rho_calc * Imcd * Km(k) * Ki(k)) / Lt;   % [24]
    Vp(k) = (rho_calc * Imcd * Ks(k) * Ki(k)) / Lt;   % [29]

    fprintf('\n--- Direcao %d (D = %.2f m, n = %d) ---\n', k, Dk, nk);
    fprintf('  kii = %.4f | Ki = %.4f\n', kii(k), Ki(k));
    fprintf('  Km  = %.4f  ->  Vt = %.2f V\n', Km(k), Vt(k));
    fprintf('  Ks  = %.4f  ->  Vp = %.2f V\n', Ks(k), Vp(k));
end

%% -----------------------------------------------------------------------
%% 7) Resistencia de aterramento (Sverak) - para referencia
%% -----------------------------------------------------------------------
A_area = 300*150;
Rg = rho_a*((1/Lt_horizontais)+(1/sqrt(20*A_area))*(1+(1/(1+H*sqrt(20/A_area)))));
fprintf('\nResistencia de aterramento (Sverak, rho_a) = %.4f ohm\n', Rg);

%% -----------------------------------------------------------------------
%% 8) Piores casos (maiores valores encontrados nas duas direcoes)
%% -----------------------------------------------------------------------
[Vt_max, idx_t] = max(Vt);
[Vp_max, idx_p] = max(Vp);

fprintf('\n=========================================================\n');
fprintf('Tensao de toque calculada (pior caso, direcao %d) = %.2f V\n', idx_t, Vt_max);
fprintf('Tensao de toque PERMISSIVEL                       = %.2f V\n', Valor_permissivel_toque);
if Vt_max <= Valor_permissivel_toque
    fprintf('  -> OK: dentro do limite\n');
else
    fprintf('  -> ATENCAO: tensao de toque calculada EXCEDE o limite permissivel!\n');
end

fprintf('\nTensao de passo calculada (pior caso, direcao %d)  = %.2f V\n', idx_p, Vp_max);
fprintf('Tensao de passo PERMISSIVEL                        = %.2f V\n', Valor_permissivel_passo);
if Vp_max <= Valor_permissivel_passo
    fprintf('  -> OK: dentro do limite\n');
else
    fprintf('  -> ATENCAO: tensao de passo calculada EXCEDE o limite permissivel!\n');
end
fprintf('=========================================================\n');
