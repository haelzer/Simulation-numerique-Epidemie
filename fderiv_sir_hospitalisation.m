function f=fderiv_sir_hospitalisation(t,x)

% voici un début de fonction dérivée, qu'il vous faut développer
% fonction dérivée pour le problème d'épidémie non spatialisée

global a b c d nlits ho hr hc
global iconf  ivariant betamax  % si iconf=1, on declenche le confinement
% si ivariant=1, on introduit le variant

alpha=1;  % coef de confinement, = 1 si pas de confinement
 
if iconf==1  % traitement du confinement
    if t>=30
       alpha=0.1;  % permier confinement sévère
    end
    if t>=60
       alpha=1; % seconde période: pas de confinement
    end
    if t>=100
       alpha=0.3; % troisième période: confinement modéré
    end
    if t>=150
       alpha=0.7; % quatrième période: confinement léger
    end


 end
   
beta=1;
tvar=200;
m=(betamax-1)/50;
r=1-m*tvar;
if ivariant==1
    if ((t>=tvar) && (t<=(tvar+50)))
        beta=m*t+r;
    end
    if t>(tvar+50)
        beta=betamax;
    end


end


%    Et puis on renseigne le vecteur derivée, à vous de le faire
nlitsdisp=max(0,nlits-x(3));
f(1)=-a*alpha*beta*x(1)*x(2);
f(2)=a*beta*alpha*x(1)*x(2)-b*x(2)-c*x(2)-min(ho*x(2),nlitsdisp);
f(3)=b*x(2)+hr*x(5);
f(4)=c*x(2)+hc*x(5);
f(5)=min(ho*x(2),nlitsdisp)-hc*x(5)-hr*x(5);
f=f';
