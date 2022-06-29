% programme principal pour le problème de l'épidémie
% SIRD avec hospitalisation


clear all
close  all
global a b c d nlits ho hr hc
global iconf  ivariant  betamax % les cas de confinement et de variant, activés si =1

% Modèle SIR classique avec décédés

iconf=0 ;   % confinement activé dans fderiv si iconf=1
ivariant=0  ;  %  variant activé dans fderiv si ivariant=1
betamax=5;  % sur-contagion due au variant, à utiliser dans fderiv


%valeurs des coefficients (influence à étudier):
sains_ini=0.995; % pourcentage de la population totale
infect_ini=0.005;
nlits=0.0005;
a=1.1e-1; % coef infection sains
b=6e-3 ; % coef guerison
c=0.02 ; % coef mortalité
d=0.; % coef resucceptibilité
ho=0.0005; % pourcentage des cas graves à hospitaliser
hr=0.1;
hc=0.1; % coef mortalite hopital
x0=[sains_ini,infect_ini,0,0,0]; % condition initiale

pastemps=1;
tmax=365;% sur1 an
intervalle_temps=[0, tmax]; 

[t,X]=Solveur(@fderiv_sir_hospitalisation,intervalle_temps,x0,pastemps,2);  % appel au solveur écrit au TP précédent

%Dessin 
figure(1)
plot(t,X,'linewidth',2);
grid on
title('Epidémie', 'FontSize', 14)
legend(' sains', 'infectés', 'rétablis', 'décédés ', 'hospitalisés', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 14)
ylabel('Populations %', 'FontSize', 14)
