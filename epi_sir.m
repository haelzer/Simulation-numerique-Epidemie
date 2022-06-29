% programme principal pour le problème de l'épidémie
% SIRD simple

% Voici le programme principal , avec l'appel à votre solveur 
% Quand vous rajouterez l'hopital, il faudra le modifier un peu 

clear all
close  all
global a b c  % les 3 coefs cinétiques
global iconf  ivariant  betamax % les cas de confinement et de variant, activés si =1

% Modèle SIR classique avec décédés

iconf=0 ;   % confinement activé dans fderiv si iconf=1
ivariant=0  ;  %  variant activé dans fderiv si ivariant=1
betamax=5;  % sur-contagion due au variant, à utiliser dans fderiv


%valeurs des coefficients (influence à étudier):
sains_ini=0.995; % pourcentage de la population totale
infect_ini=0.005;
a=1.e-1;   % coef infection sains
b=6e-3 ;   % coef guerison 
c=2e-2  ;   % coef mortalité


x0=[sains_ini,infect_ini,0, 0]; % condition initiale, 4 classes

pastemps=1;
tmax=365;% sur1 an
intervalle_temps=[0, tmax]; 

[t,X]=Solveur(@fderiv_sir,intervalle_temps,x0,pastemps,2);  % appel au solveur écrit au TP précédent

%Dessin 
figure(1)
plot(t,X,'linewidth',2);
grid on
title('Epidémie', 'FontSize', 14)
legend(' sains', 'infectés', 'rétablis', 'décédés ', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 14)
ylabel('Populations %', 'FontSize', 14)


