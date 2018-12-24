%----------DIN�MICA DE GASOS-------------

%RESOLUCI� ANAL�TICA DE PROBLEMA COMBINAT. 
%TUB AMB FLUID INTERIOR I CONVECCI� EXTERIOR.



%%-----------C�LCULS PREVIS---------------

%Discretitzaci� del tub, condicions d'entrada...
clear all;
Di = 0.01; ri = Di/2; L = 7.659; epsilon = 0.004; R = 287;

gamma = 1.4;
Min = 0.1;
Tin = 400;
vin = sqrt(gamma*R*Tin)*Min; 
pin = 5e5; 
Sup = pi*ri^2;
rhoin = pin/(287*Tin);
       
massa = Sup*vin*rhoin;

Tt = 300; %En aquest cas ens donen Ttub;

N = 50; %numero de VC



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
T_r = zeros(N,1);
Sgen = zeros(N,1);
Sgen1_v = zeros(N,1);
Sgen2_v = zeros(N,1);
det_v = zeros(N,1);
P1_v = zeros(N,1);
P2_v = zeros(N,1);
T1_v = zeros(N,1);
T2_v = zeros(N,1);
v1_v = zeros(N,1);
v2_v = zeros(N,1);
fun = @(var) 1022 - 0.166*var + 3.5025e-4*var.^2; %funcio del Cp (per Cpr)

Trs = 300; %Suposem la temperatura de recuperaci� inicial

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
y = 1;

%% Resoluci� del problema
    
for i = 1:N
    i 
    T_s(i+1) = T(i);
    P_s(i+1) = P(i);
    V_s(i+1) = v(i);
    rho_s(i+1) = rho(i);
    dif = 10e10; %Valor arbitrari per poder iterar un altre cop
    
    while dif > 1e-2
    y = 1;
    dif2 = 30;  %difer�ncia inicial entre Tr* i Tr
    
        while dif2 > 1e-2 %Iteracions per trobar el valor de Tr
 
    Ti = 0.5*(Tt+Tin)+0.22*(Trs-Tin); %Ti = temperatura de referencia, Trs anira canviant a cada iteraci�;
    vi = 0.5*(v(i)+v(i+1)); 
    Pi = 0.5*(P(i)+P(i+1));
    rhoi = Pi/(287*Ti);
    
    % C�lcul de mu i mu_w.
    
    mu = (2.5393e-5*sqrt(Ti/273.15))/(1+122/Ti);
    mu_w = (2.5393e-5*sqrt(Tt/273.15))/(1+122/Tt);
    lambda = 3.807e-3+7.4e-5*Ti;
    Cpi = 1022 - 0.166*Ti + 3.5025e-4*Ti^2;
    Cpr = integral(fun,Tin,Trs)/(Trs-Tin);
    
    %Grup adimensionals
    
    Re = rhoi*vi*Di/mu;
    
    Pr = mu*Cpi/lambda;
    
    Gz = (pi*Di)/(4*L)*Re*Pr;
    
    r = Pr^(1/3);
    
    
    %Tornem a calcular Tr i la comparem amb Trs
    Tr = Tin + r*vi^2/(2*Cpr);
    
    dif2 = abs(Tr-Trs);
    Trs = Tr;
    
    dif2(y) = dif2;
        end
    
    %Aqui ja hem calculat el valor real de Tr, i per tant de Ti.
    %Guardem les dades en vectors per comprobar resultats
    
    Cp(i) = Cpi;
    Re_v(i) = Re;
    Pr_v(i) = Pr;
    Gz_v(i) = Gz;
    mu_v(i) = mu;
    mu_w_v(i) = mu_w;
    
    
    % C�lcul del coef. de fricci�

    
    if Re < 2000
        f = 16 * Re^-1;
     elseif (5e3 < Re) && (Re < 3e4) 
        f = 0.096 * Re^-0.25;
     else
        f = 0.078 * Re^-0.2;
    end
  
    
    
    %C�lcul del Nussel i de alfa
    
    C = 0.023; m=0.8; n=0.4; K = 1;
  
    Nu = C*Re^m*Pr^n*K;
    alfa = lambda*Nu/Di;
    Nu_v(i) = Nu;
    f_v(i) = f;
    alfa_v(i) = alfa;
    
    q(i) = alfa*(Tt-Tr)*pi*Di*delta_x; %podem calcular el flux de calor amb Tr;
    
    %Resolem el sistema d'equacions
     
    syms pres;
    syms Temp;
    syms vel;
    syms dens;
    

    A_v = massa + f*rhoi*vi/4*pi*Di*delta_x;
    B_v = Sup;
    C_v = Sup*P(i)+(massa-f*rhoi*vi/4*pi*Di*delta_x)*v(i);
    
    A_t = massa*Cpi + 0.5*(alfa*Di*pi*delta_x);
    B_t = 0.5*massa + (r*alfa*Di*pi*delta_x)/(4*Cpi);
    C_t = (massa*Cpi - alfa*Di*delta_x*pi*0.5)*T(i)+(0.5*massa-(r*alfa*Di*pi*delta_x)/(4*Cpi))*v(i)^2+ alfa*Tt*Di*pi*delta_x;
    
    A = A_v*A_t*Sup - B_v*B_t*massa*R;
    B = C_v*A_t*Sup;
    C = B_v*C_t*massa*R;
    
    %Ara tenim l'equaci� quadr�tica
    syms vel;
    eqn1 = A*vel^2 - B*vel + C == 0;
    
    det = B^2-4*A*C;
    det_v(i) = det;
    
    if det < 0 %Discriminant negatiu, no te soluci� f�sica
        break;
    else
    v_sol = solve(eqn1, vel); %Guardem les dues solucions a v_sol
    end
    
    v1 = double(v_sol(1)); 
    v2 = double(v_sol(2));
    v1_v(i) = v1;
    v2_v(i) = v2;
    
    % Calculem totes les propietats amb les dues velocitats.
    
    Vol = Sup*delta_x;
    P1 = (C_v-A_v*v1)/B_v;
    T1 = (C_t-B_t*v1^2)/A_t;
    P2 = (C_v-A_v*v2)/B_v;
    T2 = (C_t-B_t*v2^2)/A_t;
    
    P1_v(i) = P1; %Guardem en vectors per poder veure-ho desrp�s si cal
    P2_v(i) = P2;
    T2_v(i) = T2;
    T2_v(i) = T1;
    
    %Trobem l'entropia generada amb les dues velocitats.
    
    Sgen1 = 1/Vol*(massa*(Cpi*log(T1/T(i))-R*log(P1/P(i))) - q(i)*Di*pi*delta_x/Tt);
    Sgen2 = 1/Vol*(massa*(Cpi*log(T2/T(i))-R*log(P2/P(i))) - q(i)*Di*pi*delta_x/Tt);
    
    Sgen1_v(i) = Sgen1; %tornem a guardar per si ac�s. 
    Sgen2_v(i) = Sgen2;
    
    if Sgen1 < 0 && Sgen2 < 0  %Si les dues s�n negatives, no t� soluci�.
        break;
    end
    
    if abs(Sgen1) > abs(Sgen2) %Ens quedem amb la que genera menys entropia
        v(i+1) = v2;
        T(i+1) = T2;
        P(i+1) = P2;
        Sgen(i) = Sgen2;
    else
        v(i+1) = v1;
        T(i+1) = T1;
        P(i+1) = P1;
        Sgen(i) = Sgen1;
    end
    
    rho(i+1) = P(i+1)/(R*T(i+1));
    
    %Guardem i comprovem les difer�ncies entre suposat-calculat.
    
    difvector(1) = abs(v(i+1)-V_s(i+1));
    difvector(2) = abs(T(i+1)-T_s(i+1));
    difvector(3) = abs(P(i+1)-P_s(i+1));
    
    dif = max(difvector);

    %Canviem el vector suposat pel calculat i tornem a iterar.
    
    V_s(i+1) = v(i+1);
    P_s(i+1) = P(i+1);
    T_s(i+1) = T(i+1);      
    x(i); %Imprimim x(i) Per saber en quin punt del conducte es torba el programa
   end
    
end
Q = sum(q)
T_out = T(N+1)
P_out = P(N+1)
v_out = v(N+1)
Re_out = Re_v(N)
Pr_out = Pr_v(N)




%% Plots per les diferents variables
figure;
plot(x,T);
title('Temperatura en funci� de x');

figure;
plot(x,P);
title('Pressi� en funci� de x');

figure;
plot(x,v);
title('Velocitat en funci� de x');

figure;
plot(x,rho);
title('Densitat en funci� de x');


