% Programme de test des m?thodes d'Euler et Euler modifi?
% C'est le programme qui sert ? tester votre fonction "solveur"
% On a commenc? ? le programmer pour vous, pour vous aider

% Essai sur l'?quation y'=2ty  sur [0,1], dont on connait bien la solution
Y0=1;       % cond initiale
a=0;
b=1;
Tspan=[a,b];
hEuler=0.05;   % Pas pour les m?thodes d'Euler
paraODE=hEuler;
options = odeset('maxstep',paraODE); % Options pour ode45, pour avoir lun pas de m?me ordre (max) 
% mais ce pas est choisi pa

% Comparaison de 3 m?thodes, la premi?re est celle de la bibliotheque
% Matlab
[x0,y0]=ode45(@derivtest1,Tspan,Y0,options);     % ode45 (Runge Kutta Matlab)
[x1,y1]=Solveur(@derivtest1,Tspan,Y0,hEuler,1);    %Euler simple
[x2,y2]=Solveur(@derivtest1,Tspan,Y0,hEuler,2);   % Euler Modifi? (RK2)
[x3,y3]=Solveur(@derivtest1,Tspan,Y0,hEuler,3) ;   % RK4

% la solution exacte
xx=a:hEuler:b;
yy=exp(xx.^2);

% pour pouvoir comparer deux y(x), il faut que les tableaux aient le m?me vecteur x: d'o? une interpolation
% voir l'aide sur cette fonction Matlab. 
y_exact0=interp1(xx,yy,x0);  % pour le x0 renvoy? par la m?thode Matlab qui n'est peut-?tre pas le m?me x que celui que l'on construit nous m?me
y_exact1=interp1(xx,yy,x1);
y_exact2=interp1(xx,yy,x2);
y_exact3=interp1(xx,yy,x3);

ermax0=max(abs(y0-y_exact0));
ermax1=max(abs(y1-y_exact1));
ermax2=max(abs(y2-y_exact2));
ermax3=max(abs(y3-y_exact3));

fprintf ( 'Voici les erreurs max des m?thodes:  \n')
fprintf ( 'erreurs max:  ODE45   %6.3f  \n', ermax0)
fprintf ( 'erreurs max:  Euler expl   %6.3f  \n', ermax1)
fprintf ( 'erreurs max:  Euler mod   %6.3f  \n', ermax2)
fprintf ( 'erreurs max:  RK4    %6.3f  \n', ermax3)

figure(1)
hold on

plot(xx,yy,'r*')
plot(x0,y0(:,1),'y')
plot(x1,y1(:,1),'b')
plot(x2,y2(:,1),'g')
plot(x3,y3(:,1),'p')

legend({'Analytique','ODE 45', 'Euler expl', 'Euler mod', 'RK4'}, 'Location','northwest','Orientation', ...
    'horizontal')
title( 'M?thodes de r?solution de y''=2ty')
grid
hold off

% Essais sur l'?quation y"= ........ a vous de jouer.. 



