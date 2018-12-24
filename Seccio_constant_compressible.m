%----------DIN�MICA DE GASOS-------------

%RESOLUCI� ANAL�TICA DE PROBLEMA COMBINAT. 
%TUB AMB FLUID INTERIOR I CONVECCI� EXTERIOR.



%%-----------C�LCULS PREVIS---------------

%Discretitzaci� del tub, condicions d'entrada...
clear all;
Di = 0.01; ri = Di/2; L = 0.07; epsilon = 0.004; R = 287;

gamma = 1.4;
Min = 0.3;
Tin = 400;
vin = sqrt(gamma*R*Tin)*Min; 
pin = 5e5; 
Sup = pi*ri^2;
rhoin = pin/(287*Tin);
       
massa = Sup*vin*rhoin;

Tt = 300; %En aquest cas ens donen Ttub;

N = 100; %numero de VC

delta_x = L/N;

x = 0:delta_x:L;

%Definicio de tots els vectors

T_s = zeros(N+1,1);
V_s = zeros(N+1,1);
P_s = zeros(N+1,1);
rho_s = zeros(N+1,1);
v = zeros(N+1,1);
T = zeros(N+1,1);
P = zeros(N+1,1);
rho = zeros(N+1,1);
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

%% Resoluci� del problema
    
for i = 1:N

    T_s(i+1) = T(i); %Suposem variables de sortida = entrada.
    P_s(i+1) = P(i);
    V_s(i+1) = v(i);
    rho_s(i+1) = rho(i);
    dif = 10e10; %Valor arbitrari per poder iterar un altre cop
    
    while dif > 1e-5

        
    Ti = 0.5*(T(i)+T(i+1)); % Funcions mitges del VC.
    vi = 0.5*(v(i)+v(i+1)); 
    Pi = 0.5*(P(i)+P(i+1));
    rhoi = Pi/(287*Ti);
    
    %C�lcul dels coeficient alfa, r i f
    
    [alfa, r, f, Cpi] = compressible(Ti, Tt, vi, Pi, Di); 
    
    %Un cop calculats seguim amb el c�lcul de Tr
    
    Tr = Tin + r*vi^2/(2*Cpi);
 
    f_v(i) = f;
    alfa_v(i) = alfa;
    Cpi = 1022 - 0.166*Ti + 3.5025e-4*Ti^2;
    q(i) = alfa*(Tt-Tr)*pi*Di*delta_x; %podem calcular el flux de calor amb Tr;
    
    %Resolem el sistema d'equacions

    A_v = massa + f*rhoi*abs(vi)*pi*Di*delta_x/4;
    B_v = Sup;
    C_v = Sup*P(i)+(massa-f*rhoi*abs(vi)*pi*Di*delta_x/4)*v(i);
    
    A_t = massa*Cpi + 0.5*alfa*Di*pi*delta_x;
    B_t = 0.5*massa + (r*alfa*Di*pi*delta_x)/(4*Cpi);
    C_t = (massa*Cpi - alfa*Di*delta_x*pi*0.5)*T(i)+(0.5*massa-(r*alfa*Di*pi*delta_x)/(4*Cpi))*v(i)^2+ alfa*Tt*Di*pi*delta_x;
    
    A = A_v*A_t*Sup - B_v*B_t*massa*R;
    B = C_v*A_t*Sup;
    C = B_v*C_t*massa*R;

    
    %Ara tenim l'equaci� quadr�tica
   
    
    det = B^2-4*A*C;
    det_v(i) = det;
    
    if det < 0 %Discriminant negatiu, no te soluci� f�sica
        error('El determinant �s negatiu');
    end
    
    v1 = (B+sqrt(B^2-4*A*C))/(2*A); 
    v2 = (B-sqrt(B^2-4*A*C))/(2*A);
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
    T1_v(i) = T1;
    
    %Trobem l'entropia generada amb les dues velocitats.
   
    Sgen1 = 1/Vol*(massa*(Cpi*log(T1/T(i))-R*log(P1/P(i))) - q(i)*Di*pi*delta_x/Tt);
    Sgen2 = 1/Vol*(massa*(Cpi*log(T2/T(i))-R*log(P2/P(i))) - q(i)*Di*pi*delta_x/Tt);
    
    Sgen1_v(i) = Sgen1; %tornem a guardar per si ac�s. 
    Sgen2_v(i) = Sgen2;
    
    if Sgen1 < 0 && Sgen2 < 0  %Si les dues s�n negatives, no t� soluci�.
        error('Les dues entropies s�n negatives');
    end
    
    if isreal(Sgen1) == 0 && isreal(Sgen2) == 0 %Fem que si hi ha una netropia im. digui quina �s.
        error('Les dues entropies s�n imaginaries');
    end
    
    if (isreal(Sgen1) == 0 && Sgen2 < 0) || (isreal(Sgen2) == 0 && Sgen1 < 0)
        disp('Entropies no es poden resoldre');
    end
        
    if isreal(Sgen1) == 0 && Sgen2 > 0 %Comprovem si Sgen1 es real o no
        v(i+1) = v2;
        T(i+1) = T2;
        P(i+1) = P2;
        Sgen(i) = Sgen2;
        disp('Entropia 1 �s imaginaria, entropia 2 positiva');
    end
    if isreal(Sgen2) == 0 && Sgen1 > 0 %Comprovem si Sgen2 es real o no
        v(i+1) = v1;
        T(i+1) = T1;
        P(i+1) = P1;
        Sgen(i) = Sgen1;
        disp('Entropia 2 �s imaginaria, entropia 1 positiva');
    end
    if abs(Sgen1) > abs(Sgen2) %Cas en que les dues s�n positives
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
Q = sum(q);
T_out = T(N+1);
P_out = P(N+1);
v_out = v(N+1);
Re_out = Re_v(N);
Pr_out = Pr_v(N);



%% Verificaci� del codi (Guillem)

%4.1 Conservaci� de la massa

m_out = v(N+1)*rho(N+1)*Sup;
Err_mass = abs(massa-m_out)*100/massa;
if Err_mass<0.1
    disp('Conservaci� massa OKEY')
end

% 4.2 Conservaci� momentum
Momentum_in = massa*vin+pin*Sup; %Momentum entrada

%C�lcul fregament
Fregament = 0;
for i=1:N
    vi = (v(i)+v(i+1))/2;
    rho_i = (rho(i)+rho(i+1))/2;
    
    mu_i = 2.5393*10^(-5) * sqrt(T(i)/273.15) / (1+(122/T(i))); %Viscositat din�mica
    
    Fregament = Fregament + f_v(i)*rho_i*vi^2*pi*Di*delta_x/2;   
end

Momentum_out = massa*v(N+1)+P(N+1)*Sup;
Err_mom = abs(Momentum_in-Fregament-Momentum_out)*100/Momentum_in;
if Err_mom<0.1
    disp('Conservaci� momentum OKEY')
end


% 4.3 Conservaci� energia

%Energia entrada
Cp_in = 1022-0.1626*Tin+3.5025*10^(-4)*Tin^2;
E_in = massa*(Cp_in*Tin+(vin^2)/2); 

%Energia sortida
Cp_out = 1022-0.1626*T(N+1)+3.5025*10^(-4)*T(N+1)^2;
E_out = massa*(Cp_out*T(N+1)+(v(N+1)^2)/2); 
%C�lcul del calor
Q = 0;
 
for i=1:N
        V_i = (v(i)+v(i+1))/2;
        T_i = (T(i)+T(i+1))/2;
        P_i = (P(i)+P(i+1))/2;

        Cp_i = 1022-0.1626*T_i+3.5025*10^(-4)*T_i^2;

        [alfa, r, f, Cpi] = compressible(T_i, Tt, V_i, P_i, Di);
        Tr = T_i + r*V_i^2/(Cp_i*2); %Temperatura de recuperaci�
        Q = Q + alfa*(Tt-Tr)*pi*Di*delta_x;   
end


Err_ene = abs(E_in-E_out+Q)*100/E_in;
if Err_ene<0.1
    disp('Conservaci� energia OKEY')
end

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


