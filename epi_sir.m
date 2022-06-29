% programme principal pour le probl�me de l'�pid�mie
% SIRD simple

% Voici le programme principal , avec l'appel � votre solveur 
% Quand vous rajouterez l'hopital, il faudra le modifier un peu 

clear all
close  all
global a b c  % les 3 coefs cin�tiques
global iconf  ivariant  betamax % les cas de confinement et de variant, activ�s si =1

% Mod�le SIR classique avec d�c�d�s

iconf=0 ;   % confinement activ� dans fderiv si iconf=1
ivariant=0  ;  %  variant activ� dans fderiv si ivariant=1
betamax=5;  % sur-contagion due au variant, � utiliser dans fderiv


%valeurs des coefficients (influence � �tudier):
sains_ini=0.995; % pourcentage de la population totale
infect_ini=0.005;
a=1.e-1;   % coef infection sains
b=6e-3 ;   % coef guerison 
c=2e-2  ;   % coef mortalit�


x0=[sains_ini,infect_ini,0, 0]; % condition initiale, 4 classes

pastemps=1;
tmax=365;% sur1 an
intervalle_temps=[0, tmax]; 

[t,X]=Solveur(@fderiv_sir,intervalle_temps,x0,pastemps,2);  % appel au solveur �crit au TP pr�c�dent

%Dessin 
figure(1)
plot(t,X,'linewidth',2);
grid on
title('Epid�mie', 'FontSize', 14)
legend(' sains', 'infect�s', 'r�tablis', 'd�c�d�s ', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 14)
ylabel('Populations %', 'FontSize', 14)


