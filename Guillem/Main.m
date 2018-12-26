%Dinàmica de Gasos i Transferència de Calor i Massa
%Projecte voluntari - Tub amb flux compressible
%2018-19 Q1

close all;
clearvars; 

%% 1. Entrada dades -----------------------------------------------------

%Dades geomètriques
L = 0.05; %Longitud tub [m]
R = 0.005; %Radi tub [m]
epsilon = 0.004; %Relative roughness [m]

%Dades numèriques

% for n=1:10
Nvc = 1000; %Volums de control
% N(n) = Nvc;
Delta = 1e-10; %Delta per convergencia



%Dades flux d'entrada
P_in = 5e5; %Pressió d'entrada [Pa]
T_in = 400; %Temperatura d'entrada [K]

% for m=1:7
M_in = 0.6; %Mach d'entrada
% M(m)=M_in;

%Dades físiques
R_gas = 287; %Cte gas ideal [Pa*m^3/(kg*K)]
T_wall = 300; %Temperatura tub [K]

Adiabatic = false; %Si fem el cas adiabatic li donem valor true
%% 2. Càlculs previs ----------------------------------------------------

Delta_z = L/Nvc; %Longitud de cada volum de control (uniformes)
z = 0:Delta_z:L; %Posició de cada node [m]
rho_in = P_in/(R_gas*T_in); %Densitat d'entrada [kg/m^3]
V_in = M_in*sqrt(1.4*R_gas*T_in); %Velocitat d'entrada [m/s]
      
S = pi*R^2; %Secció del tub [m^2]
Per = pi*2*R; %Perímetre tub [m]
m_in = rho_in*V_in*S; %Cabal màssic d'entrada [kg/s]

%Creem vectors per les nostres variables a cada node:
V = zeros(1, Nvc+1) + V_in;
P = zeros(1, Nvc+1) + P_in;
Rho = zeros(1, Nvc+1) + rho_in;
T = zeros(1, Nvc+1) + T_in;
Mach = zeros(1, Nvc+1) + M_in;
Sgen = zeros(1, Nvc+1);
s = zeros(1, Nvc+1);


%% 3. Resolució step by step --------------------------------------------

for i=1:Nvc %Visitarem tots els VdC
    
    %Al proper node li donem el valor inicial que té el node actual
    V(i+1) = V(i);
    T(i+1) = T(i);
    Rho(i+1) = Rho(i);
    P (i+1) = P(i);
    
    for j=1:1000 %Màxim de 1000 iteracions a cada VdC

        %Càlculs valors mitjans al VdC
        V_i = (V(i)+V(i+1))/2;
        T_i = (T(i)+T(i+1))/2;
        Rho_i = (Rho(i)+Rho(i+1))/2;
        P_i =(P(i)+P(i+1))/2;
        
        % Cas adiabàtic
        alpha = 0; 
        r = 0;
        mu_i = (1.458*10^(-6)*T_i^(1.5))/(T_i+110.4); %Viscositat dinàmica
        Cp_i = 1034.09-2.849*10^(-1)*T_i+7.817*10^(-4)*T_i^2-4.971*10^(-7)*T_i^3+1.088*10^(-10)*T_i^4;
       
        % Cas NO adiabàtic
        if Adiabatic~=1 
            [r, alpha, mu_i, Cp_i] = Recuperacio(T_i, T_wall, P_i, V_i, R, Delta); %Coeficients
        end      
        
        Tr = T_i + r*V_i^2/(Cp_i*2); %Temperatura de recuperació
        
        f_i = Rugositat (V_i, R, epsilon, Rho_i, mu_i); %Skin friction coefficient
        %Cp_i = 1022-0.1626*T_i+3.5025*10^(-4)*T_i^2;
        
        
                
        %Eq. del momentum
        %Tenim Av*v(i+1)+Bv*p(i+1)-Cv = 0
        Av = m_in + f_i*Rho_i*abs(V_i)*Per*Delta_z/4;
        Bv = S;
        Cv = S*P(i)+(m_in - f_i*Rho_i*abs(V_i)*Per*Delta_z/4)*V(i);

        %Eq. de l'energia
        %At*T(i+1)+Bt*v^2(i+1)-Ct = 0
        At = m_in*Cp_i + alpha*Per*Delta_z/2;
        Bt = m_in/2 + r*alpha*Per*Delta_z/(4*Cp_i);
        Ct = (m_in*Cp_i - alpha*Per*Delta_z/2)*T(i)  +  (m_in/2 - r*alpha*Per*Delta_z/(4*Cp_i))*V(i)^2  +  alpha*T_wall*Per*Delta_z; 

        %Introduint conservació de massa i eq d'estat en l'eq del momentum, i
        %després subsituint T(i+1) a l'eq de l'energia, podem reagrupar:
        %A*v(i+1)^2-B*v(i+1)+C = 0
        A = Av*At*S - Bv*Bt*m_in*R_gas;
        B = Cv*At*S;
        C = Bv*Ct*m_in*R_gas;

        %Resolem l'equació per trobar v(i+1)
        %S el discriminant és negatiu no hi ha solució física
        if B^2-4*A*C < 0
            error = i
            break
        end

        %Si és positiu, tenim dues solucions
        v_1 = (B+sqrt(B^2-4*A*C))/(2*A);
        v_2 = (B-sqrt(B^2-4*A*C))/(2*A);

        %Evaluem l'entropia generada d'aquestes
        p_1 = (Cv-Av*v_1)/Bv;
        T_1 = (Ct-Bt*v_1^2)/At;
        S_gen_1 = (1/(S*Delta_z))*(m_in*(Cp_i*log(T_1/T(i))-R_gas*log(p_1/P(i))) - alpha*(T_wall-Tr)*Per*Delta_z/T_wall);

        p_2 = (Cv-Av*v_2)/Bv;
        T_2 = (Ct-Bt*v_2^2)/At;
        S_gen_2 = (1/(S*Delta_z))*(m_in*(Cp_i*log(T_2/T(i))-R_gas*log(p_2/P(i))) - alpha*(T_wall-Tr)*Per*Delta_z/T_wall);

        %Ja estem en condicions de triar quina és la solució amb sentit físic

        if (S_gen_1>=0 && S_gen_2<0) || imag(S_gen_2)~=0 %Ens quedem amb la 1
            v_calc = v_1;
            P_calc = p_1;
            T_calc = T_1;
            Rho_calc = P_calc/(R_gas*T_calc);
            Sgen(i+1) = S_gen_1;
        end

        if (S_gen_2>=0 && S_gen_1<0) || imag(S_gen_1)~=0 %Ens quedem amb la 2
            V_calc = v_2;
            P_calc = p_2;
            T_calc = T_2;
            Rho_calc = P_calc/(R_gas*T_calc);
            Sgen(i+1) = S_gen_2;
        end
        
        if S_gen_1>=0 && S_gen_2>=0  %Si les dues tenen entropia positiva amb quina ens quedem?
            
            %Ens quedem amb la ''mes propera'', de moment no ha d'haver-hi ones
            %de xoc
            if abs(p_2-P_i) < abs(p_1-P_i)
                 V_calc = v_2;
                 P_calc = p_2;
                 T_calc = T_2;
                 Rho_calc = P_calc/(R_gas*T_calc); 
                 Sgen(i+1) = S_gen_2;
            end
            
            if abs(p_2-P_i) > abs(p_1-P_i)
                 V_calc = v_1;
                 P_calc = p_1;
                 T_calc = T_1;
                 Rho_calc = P_calc/(R_gas*T_calc);
                 Sgen(i+1) = S_gen_1;
            end
        end
     
        %Comprovem convergència
        if abs(V_calc-V(i+1))<Delta && abs(T_calc-T(i+1))<Delta && abs(P_calc-P(i+1))<Delta && abs(Rho_calc-Rho(i+1))<Delta
            V(i+1) = V_calc;
            P(i+1) = P_calc;
            T(i+1) = T_calc;
            Rho(i+1) = Rho_calc;
            c = sqrt(1.4*R_gas*T(i+1));
            Mach(i+1) = V(i+1)/c;
            break
            %Si ha complert convergència, que passi al següent node
        end
        
        %Si no, tenim nous valors per la nova iteració
        V(i+1) = V_calc;
        P(i+1) = P_calc;
        T(i+1) = T_calc;
        Rho(i+1) = Rho_calc;
        c = sqrt(1.4*R_gas*T(i+1));
        Mach(i+1) = V(i+1)/c;
        s(i+1) = (alpha*(T_wall-Tr)*Per*Delta_z/T_wall + Sgen(i+1)*Delta_z*S)*(1/m_in)+s(i);
    end
        
    
end

%% 4. Verificació del codi -----------------------------------------------
% Verificarem balanços globals
n=1; m=1;
% 4.1 Conservació de la massa
m_out = V(Nvc+1)*Rho(Nvc+1)*S;
Err_mass(n,m) = abs(m_in-m_out)*100/m_in;
if Err_mass(n,m)<0.1
    disp('Conservació massa OKEY')
end


% 4.2 Conservació momentum
Momentum_in = m_in*V_in+P_in*S; %Momentum entrada

%Càlcul fregament
Fregament = 0;
for i=1:Nvc
    V_i = (V(i)+V(i+1))/2;
    Rho_i = (Rho(i)+Rho(i+1))/2;
    T_i = (T(i)+T(i+1))/2;
    P_i = (P(i)+P(i+1))/2;
    
    [r,alpha,mu_i,Cp] = Recuperacio (T_i, T_wall, P_i, V_i, R, Delta);
    f_i = Rugositat(V_i, R, epsilon, Rho_i, mu_i); %Skin friction coefficient
    
    Fregament = Fregament + f_i*Rho_i*V_i^2*Per*Delta_z/2;   
end

Momentum_out = m_in*V(Nvc+1)+P(Nvc+1)*S;
Err_mom(n,m) = abs(Momentum_in-Fregament-Momentum_out)*100/Momentum_in;
if Err_mom(n,m)<0.1
    disp('Conservació momentum OKEY')
end

% 4.3 Conservació energia
%Energia entrada
[~,~,~,Cp_in] = Recuperacio(T_in,T_wall,P_in,V_in,R,Delta);
E_in = m_in*(Cp_in*T_in+(V_in^2)/2); 

%Energia sortida
[~,~,~,Cp_out] = Recuperacio(T(Nvc+1),T_wall,P(Nvc+1),V(Nvc+1),R,Delta);
E_out = m_in*(Cp_out*T(Nvc+1)+(V(Nvc+1)^2)/2); 

%Càlcul del calor
Q = 0;
if Adiabatic~=1 
    for i=1:Nvc
        V_i = (V(i)+V(i+1))/2;
        T_i = (T(i)+T(i+1))/2;
        P_i = (P(i)+P(i+1))/2;

        [r, alpha, mu, Cp_i] = Recuperacio (T_i, T_wall, P_i, V_i, R, Delta); %Coeficients
        Tr = T_i + r*V_i^2/(Cp_i*2); %Temperatura de recuperació
        Q = Q + alpha*(T_wall-Tr)*Per*Delta_z;   
    end
end

Err_ene(n,m) = abs(E_in-E_out+Q)*100/E_in;
if Err_ene(n,m)<0.1
    disp('Conservació energia OKEY')
end

% end
% end

%% 5. Impressió de resultats ---------------------------------------------

figure()
subplot(2,2,1)
plot(z,V)
xlabel ('x [m]')
ylabel('v [m/s]')
title('Velocitat')

subplot(2,2,2)
plot(z,P)
title('Pressió')
xlabel ('x [m]')
ylabel('P [Pa]')

subplot(2,2,3)
plot(z,T)
title('Temperatura')
xlabel ('x [m]')
ylabel('T [K]')

subplot(2,2,4)
plot(z,Rho)
title('Densitat')
xlabel ('x [m]')
ylabel('Rho [kg/m^3]')





figure()
subplot(2,2,1)
plot(z,Mach)
title('Mach')
xlabel ('x [m]')
ylabel('Ma')


subplot(2,2,2)
plot(z,Sgen)
title('Entropia generada')
xlabel ('x [m]')
ylabel('Sgen [w/m^3 K]')

subplot(2,2,3)
plot(z,s)
title('Entropia acumulada')
xlabel ('x [m]')
ylabel('s [J/kgK]')


% figure()
% plot(N,Err_mom)
% grid on
% xlabel('Número de volums de control')
% ylabel('Error momentum')
% legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')
% 
% figure()
% plot(N,Err_mass)
% grid on
% xlabel('Número de volums de control')
% ylabel('Error massa')
% legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')
% 
% figure()
% plot(N,Err_ene)
% grid on
% xlabel('Número de volums de control')
% ylabel('Error energia')
% legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')




