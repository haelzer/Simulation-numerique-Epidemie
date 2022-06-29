function f=fderiv_sir(t,x)

% voici un d�but de fonction d�riv�e, qu'il vous faut d�velopper
% fonction d�riv�e pour le probl�me d'�pid�mie non spatialis�e

global a b c  
global iconf  ivariant betamax  % si iconf=1, on declenche le confinement
% si ivariant=1, on introduit le variant

alpha=1;  % coef de confinement, = 1 si pas de confinement
 
if iconf==1  % traitement du confinement
    if t>=30
       alpha=0.1;  % permier confinement s�v�re
    end
    if t>=60
       alpha=1; % seconde p�riode: pas de confinement
    end
    if t>=100
       alpha=0.3; % troisi�me p�riode: confinement mod�r�
    end
    if t>=150
       alpha=0.7; % quatri�me p�riode: confinement l�ger
    end


 end
   
beta=1;
tvar=100;
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


%    Et puis on renseigne le vecteur deriv�e, � vous de le faire
% f(1)=........
% ................
f(1)=-a*alpha*beta*x(1)*x(2);
f(2)=a*beta*alpha*x(1)*x(2)-b*x(2)-c*x(2);
f(3)=b*x(2);
f(4)=c*x(2);

f=f'; % renvoie vecteur deriv�e colonne

