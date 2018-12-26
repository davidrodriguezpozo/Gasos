%----------DINÀMICA DE GASOS-------------

%RESOLUCIÓ ANALÍTICA DE PROBLEMA COMBINAT. 
%TUB AMB FLUID INTERIOR I CONVECCIÓ EXTERIOR.



%%-----------CÀLCULS PREVIS---------------

%Discretització del tub, condicions d'entrada...
clear all;
Di = 0.02; ri = Di/2; L = 5; epsilon = 0.0001;

water = true; air = false; oil = false; %Tipus de liquid

vin = 1; 
pin = 2e5; 
Tin = 20+273.15; 
Sup = pi*ri^2;


if water
rhoin = 847.2 + 1.298*Tin - 2.657e-3*Tin^2;
elseif oil
rhoin = 1164.45 - 0.4389*Tin-3.21e-4*Tin^2;
else
rhoin = pin/(287*Tin);
end
       
massa = Sup*vin*rhoin;

Tt = 95+273.15; %En aquest cas ens donen Ttub;

N = 100000; %numero de VC



delta_x = L/N;

x = 0:delta_x:L;

%Definicio de tots els vectors

S = zeros(N+1,1);
Cp = zeros(N+1,1);
A = zeros(N+1,1);
T_s = zeros(N+1,1);
V_s = zeros(N+1,1);
P_s = zeros(N+1,1);
rho_s = zeros(N+1,1);
v = zeros(N+1,1);
T = zeros(N+1,1);
P = zeros(N+1,1);
rho = zeros(N+1,1);
mu_v = zeros(N+1,1);
lambda = zeros(N+1,1);
Re_v = zeros(N,1);
Pr_v = zeros(N,1);
Nu_v = zeros(N,1);
alfa_v = zeros(N,1);
mu_w_v = zeros(N,1);
Gz_v = zeros(N,1);
f_v = zeros(N,1);
q = zeros(N,1);
%Discretitzem a tots els nodes (variables que no canvien), i omplim els
%vectors amb les variables suposades.
for i = 1:N+1
    S(i) = pi*ri^2;
    A(i) = Di*pi*delta_x;
    T_s(i) = Tin;
    V_s(i) = vin;
    P_s(i) = pin;
    rho_s(i) = rhoin;
end

for i = 1:N+1
    S(i) = pi*ri^2;
    A(i) = Di*pi*delta_x;
    T(i) = Tin;
    v(i) = vin;
    P(i) = pin;
    rho(i) = rhoin;
end


%% Resolució del problema
    
for i = 1:N
    T_s(i+1) = T(i);
    P_s(i+1) = P(i);
    V_s(i+1) = v(i);
    rho_s(i+1) = rho(i);
    dif = 10e10; %Valor arbitrari per poder iterar un altre cop
    
    while dif > 1e-10
        
    Ti = 0.5*(T(i)+T(i+1));
    Pi = 0.5*(P(i)+P(i+1));
    
    if water
    rhoi = 847.2 + 1.298*Ti - 2.657e-3*Ti^2;
    elseif oil
    rhoi = 1164.45 - 0.4389*Ti-3.21e-4*Ti^2;
    else
    rhoi = Pi/(287*Ti);
    end
    
    vi = 0.5*(v(i)+v(i+1));
    
    
    if water %Tractem el cas de l'aigua
    
       
    if Ti < 353 %Càlcul de mu
       mu = 0.9149 -1.2563e-2*Ti+6.9182e-5*Ti^2-1.9067e-7*Ti^3+2.6275e-10*Ti^4-1.4474e-13*Ti^5;
    else
       mu = 3.7471e-2-3.5636e-4*Ti + 1.3725e-6*Ti^2-2.6566e-9*Ti^3+2.5766e-12*Ti^4-1e-15*Ti^5;
    end
    
    if Tt < 353 %càlcul de mu_w
        mu_w = 0.9149 - 1.2563e-2*Tt+ 6.9182e-5*Tt^2 - 1.9067e-7*Tt^3 + 2.476e-10*Tt^41.4474e-13*Tt^5;
    else
        mu_w = 3.7471e-2 - 3.5636e-4*Tt + 1.3725e-6*Tt^2 - 2.6566e-9*Tt^3 + 2.5766e-12*Tt^4-1e-15*Tt^5;
    end
    
    lambda = -1.176 + 7.915e-3*Ti + 1.486e-5*Ti^2 - 1.317e-7*Ti^3 + 2.476e-10*Ti^4 - 1.556e-13*Ti^5;
    
    Cpi = 5648.8 - 9.140*Ti + 14.21e-3*Ti^2;
    
    elseif air %Ara tractem el cas de l'aire
    
    if Ti < 1500 %Càlcul de mu
       mu = (1.458e-6*Ti^1.5)/(Ti+110.4);
    else
       mu = (2.5393e-5*sqrt(Ti/273.15))/(1+122/Ti);
    end
   
    if Tt < 1500  %càlcul de mu_w
       mu_w = (1.458e-6*Tt^1.5)/(Tt+110.4);
    else
       mu_w = (2.5393e-5*sqrt(Tt/273.15))/(1+122/Tt);
    end
    
    lambda = (2.648e-3*sqrt(Ti))/(1+(245.4/Ti)*10^(-12/Ti));
    
    Cpi = 1034.09-2.849e-1*Ti+7.817e-4*Ti^2-4.971e-7*Ti^3+1.077e-10*Ti^4;
    
    else %Tractem el cas de l'oli
    
    mu = rhoi*exp(-16.096+(586.38/(Ti-210.65))); %càlcul de mu
   
    mu_w = rhoi*exp(-16.096+(586.38/(Tt-210.65))); %calcul de mu_w
    
    lambda = 0.116 + 4.9e-5*Ti - 1.5e-7*Ti^2;
    
    Cpi = 658 + 2.82*Ti + 8.97e-4*Ti^2;    
        
    end
    Re = rhoi*vi*Di/mu;
    
    Pr = mu*Cpi/lambda;
    
    Gz = (pi*Di)/(4*L)*Re*Pr;
    
    %Guardem les dades en vectors per comprobar resultats
    
    Cp(i) = Cpi;
    Re_v(i) = Re;
    Pr_v(i) = Pr;
    Gz_v(i) = Gz;
    mu_v(i) = mu;
    mu_w_v(i) = mu_w;
    
    
    % Càlcul del coef. de fricció
    
    if Re < 2000
        f = 16 * Re^-1;
    elseif (5e3 < Re) && (Re < 3e4) 
        f = 0.079 * Re^-0.25;
    else
        f = 0.046 * Re^-0.2;
    end
    
    
    %Càlcul del Nussel i de alfa
    
    if Re < 2000 && Gz> 10 %estem a flux laminar, tub curt. 
       C = 1.86; m = 1/3; n = 1/3; K = (Di/L)^(1/3)*(mu/mu_w)^0.14;
    
    elseif Re < 2000 && Gz < 10 %Flux Laminar, tub llarg.
       C = 3.66; m = 0; n = 0; K = 1;
    else %Flux turbulent.
        C = 0.027; m= 0.8; n = 0.33; K = (mu/mu_w)^0.14;
    end
    
    if air && Re > 2000
        C = 0.023; m=0.8; n=0.4; K = 1;
    end
    
    Nu = C*Re^m*Pr^n*K;
    alfa = lambda*Nu/Di;
    Nu_v(i) = Nu;
    f_v(i) = f;
    alfa_v(i) = alfa;
    
    %Resolem el sistema d'equacions
    
    P(i+1) = 1/Sup*(-massa*(v(i+1)-v(i))-0.5*f*rhoi*vi^2*pi*Di*delta_x+P(i)*Sup);
    T(i+1) = 1/(massa*Cpi+alfa*0.5*pi*Di*delta_x)*(alfa*(Tt-0.5*T(i))*pi*Di*delta_x+massa*Cpi*T(i)-massa*0.5*(v(i+1)^2-v(i)^2));
    %càlcul de les diferents densitats per separat.
    
    if water
    rho(i+1) = 847.2 + 1.298*T(i+1) - 2.657e-3*T(i+1)^2;
    elseif oil
    rho(i+1) = 1164.45 - 0.4389*T(i+1)-3.21e-4*T(i+1)^2;
    else
    rho(i+1) = P(i+1)/(287*T(i+1));    
    end
    
    v(i+1) = massa/(rho(i+1)*Sup);
    
    q(i) = alfa*(Tt-0.5*(T(i)+T(i+1)))*pi*Di*delta_x;
   
    

    %Guardem i comprovem les diferències entre suposat-calculat.
    
    difvector(1) = abs(v(i+1)-V_s(i+1));
    difvector(2) = abs(T(i+1)-T_s(i+1));
    difvector(3) = abs(P(i+1)-P_s(i+1));
    
    dif = max(difvector);

    %Canviem el vector suposat pel calculat i tornem a iterar.
    
    V_s(i+1) = v(i+1);
    P_s(i+1) = P(i+1);
    T_s(i+1) = T(i+1);      
    
   end
    
end
Q = sum(q);
T_out = T(N+1)-273.15;
P_out = P(N+1);
v_out = v(N+1);
Re_out = Re_v(N);
Pr_out = Pr_v(N);




%% Plots per les diferents variables
figure;
plot(x,T);
title('Temperatura en funció de x');

figure;
plot(x,P);
title('Pressió en funció de x');

figure;
plot(x,v);
title('Velocitat en funció de x');

figure;
plot(x,rho);
title('Densitat en funció de x');


