%Dinàmica de Gasos i Transferència de Calor i Massa
%Projecte voluntari - Tub amb flux compressible
%2018-19 Q1
%Guillem Vergés i Plaza

close all;
%clearvars -except Rho T V  %%% PER LESTUDI DE DELTA
clearvars

% Rho_error = Rho(2001); %%% PER LESTUDI DE DELTA
% T_error = T(2001);%%% PER LESTUDI DE DELTA
% V_error =V(2001);%%% PER LESTUDI DE DELTA

%% 1. Entrada dades -----------------------------------------------------

%Dades geomètriques
L = 0.5; %Longitud tub [m]
epsilon = 0.0001; %Relative roughness [m]

%Dades numèriques

%for n=1:100
Nvc = 2000; %Volums de control
%N = N(n);

%Estudiem la influència del factor de convergència
% Delt = linspace(1,20,120); %%% PER LESTUDI DE DELTA
% for n=1:length(Delt) %%% PER LESTUDI DE DELTA
% Delta = 10^-Delt(n); %Delta per convergencia %%% PER LESTUDI DE DELTA
% D(n) = Delta; %%% PER LESTUDI DE DELTA
% tic %%% PER LESTUDI DE DELTA

Delta = 10^(-9); %Delta per convergencia


%Dades flux d'entrada
P_in = 5e5; %Pressió d'entrada [Pa]
T_in = 600; %Temperatura d'entrada [K]

%Quants Machs d'entrada estudiarem?
Machs = 1; %%% EN AQUESTA VERSIÓ HO VALORAREM COM ONES DE XOC

for m=1:Machs
M_in = 0.049597; %Mach d'entrada per cada iteració, en provarem diferents

M = [2*L 0.7*L 0.9*L 15*L]; %%% ARA SERÀ EL VECTOR DE POSICIÓ ONA DE XOC

%Dades físiques
R_gas = 287; %Cte gas ideal [Pa*m^3/(kg*K)]
T_wall = 300; %Temperatura tub [K]

Adiabatic = true; %Si fem el cas adiabatic li donem valor true
%% 2. Càlculs previs ----------------------------------------------------

Delta_z = L/Nvc; %Longitud de cada volum de control (uniformes)
z = 0:Delta_z:L; %Posició de cada node [m]
R = z + 0.015 +1./(100*z+1.5); %Radi de cada node [m]
rho_in = P_in/(R_gas*T_in); %Densitat d'entrada [kg/m^3]
V_in = M_in*sqrt(1.4*R_gas*T_in); %Velocitat d'entrada [m/s]
      
S = pi*R.^2; %Secció del tub [m^2]
Per = pi*2.*R; %Perímetre tub [m]
for i=1:Nvc
    Theta(i) = atan((R(i+1)-R(i))/Delta_z); %Angle a cada VC
    A_i(i) = 2*pi*(R(i)+R(i+1))*0.5*Delta_z/cos(Theta(i)); %Area normal a cada VC
end
    
m_in = rho_in*V_in*S(1); %Cabal màssic d'entrada [kg/s]

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
        
        %Mirem el valor que té a cada iteració per estudiar la convergència
        V_j(j) = V_i;
        T_j(j) = T_i;
        Rho_j(j) = Rho_i;
        P_j (j) = P_i;
        
        % Cas adiabàtic
        alpha = 0; 
        r = 0;
        mu_i = (1.458*10^(-6)*T_i^(1.5))/(T_i+110.4); %Viscositat dinàmica
        Cp_i = 1034.09-2.849*10^(-1)*T_i+7.817*10^(-4)*T_i^2-4.971*10^(-7)*T_i^3+1.088*10^(-10)*T_i^4;
       
        % Cas NO adiabàtic
        if Adiabatic~=1 
            [r, alpha, mu_i, Cp_i] = Recuperacio(T_i, T_wall, P_i, V_i, R(i), Delta); %Coeficients
        end      
        
        Tr = T_i + r*V_i^2/(Cp_i*2); %Temperatura de recuperació
        
        f_i = Rugositat (V_i, R(i), epsilon, Rho_i, mu_i); %Skin friction coefficient
        %Cp_i = 1022-0.1626*T_i+3.5025*10^(-4)*T_i^2;
        
        
                
        %Eq. del momentum
        %Tenim Av*v(i+1)+Bv*p(i+1)-Cv = 0
        Av = m_in + f_i*Rho_i*abs(V_i)*A_i(i)*cos(Theta(i))/4;
        Bv = S(i+1)-A_i(i)*sin(Theta(i))/2;
        Cv = (S(i)+ A_i(i)*sin(Theta(i))*0.5)*P(i) + (m_in - f_i*Rho_i*abs(V_i)*A_i(i)*cos(Theta(i))/4)*V(i);

        %Eq. de l'energia
        %At*T(i+1)+Bt*v^2(i+1)-Ct = 0
        At = m_in*Cp_i + alpha*A_i(i)/2;
        Bt = m_in/2 + r*alpha*A_i(i)/(4*Cp_i);
        Ct = (m_in*Cp_i - alpha*A_i(i)/2)*T(i)  +  (m_in/2 - r*alpha*A_i(i)/(4*Cp_i))*V(i)^2  +  alpha*T_wall*A_i(i); 

        %Introduint conservació de massa i eq d'estat en l'eq del momentum, i
        %després subsituint T(i+1) a l'eq de l'energia, podem reagrupar:
        %A*v(i+1)^2-B*v(i+1)+C = 0
        A = Av*At*S(i) - Bv*Bt*m_in*R_gas;
        B = Cv*At*S(i);
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
        S_gen_1 = (1/(S(i)*Delta_z))*(m_in*(Cp_i*log(T_1/T(i))-R_gas*log(p_1/P(i))) - alpha*(T_wall-Tr)*Per(i)*Delta_z/T_wall);

        p_2 = (Cv-Av*v_2)/Bv;
        T_2 = (Ct-Bt*v_2^2)/At;
        S_gen_2 = (1/(S(i)*Delta_z))*(m_in*(Cp_i*log(T_2/T(i))-R_gas*log(p_2/P(i))) - alpha*(T_wall-Tr)*Per(i)*Delta_z/T_wall);

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
        
        if S_gen_1>=0 && S_gen_2>=0 && imag(S_gen_2)==0 && imag(S_gen_1)==0 %Si les dues tenen entropia positiva amb quina ens quedem?
            
            %Ens quedem amb la ''mes propera'', de moment no ha d'haver-hi ones
            %de xoc
            
            %%% ARA IMPLEMENTEM LA ONA DE XOC
            
            if M(m)==(z(i)-L/Nvc) %Si toca que hi hagi ona de xoc

                if abs(p_2-P_i) > abs(p_1-P_i)
                     V_calc = v_2;
                     P_calc = p_2;
                     T_calc = T_2;
                     Rho_calc = P_calc/(R_gas*T_calc); 
                     Sgen(i+1) = S_gen_2;
                end

                if abs(p_2-P_i) < abs(p_1-P_i)
                     V_calc = v_1;
                     P_calc = p_1;
                     T_calc = T_1;
                     Rho_calc = P_calc/(R_gas*T_calc);
                     Sgen(i+1) = S_gen_1;
                end
            end
            
            if M(m)~=(z(i)-L/Nvc)
            % Si no toca ona de xoc, ens quedem amb les variacions petites
                if v_2>v_1
                     V_calc = v_2;
                     P_calc = p_2;
                     T_calc = T_2;
                     Rho_calc = P_calc/(R_gas*T_calc); 
                     Sgen(i+1) = S_gen_2;
                end

                if v_2<v_1
                     V_calc = v_1;
                     P_calc = p_1;
                     T_calc = T_1;
                     Rho_calc = P_calc/(R_gas*T_calc);
                     Sgen(i+1) = S_gen_1;
                end
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
            
            J(i)=j;
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
        s(i+1) = (alpha*(T_wall-Tr)*A_i(i)/T_wall + Sgen(i+1)*Delta_z*S(i))*(1/m_in)+s(i);
    end
        
    
end

%% 4. Verificació del codi -----------------------------------------------
% Verificarem balanços globals

%Si només estem analitzant un cas n=m=1
n=1; 
%m=1;


% 4.1 Conservació de la massa
m_out = V(Nvc+1)*Rho(Nvc+1)*S(Nvc+1);
Err_mass(n,m) = abs(m_in-m_out)*100/m_in;
if Err_mass(n,m)<0.1
    disp('Conservació massa OKEY')
end


% 4.2 Conservació momentum
Momentum_in = m_in*V_in+P_in*S(1); %Momentum entrada

%Càlcul fregament
Fregament = 0;
for i=1:Nvc
    V_i = (V(i)+V(i+1))/2;
    Rho_i = (Rho(i)+Rho(i+1))/2;
    T_i = (T(i)+T(i+1))/2;
    P_i = (P(i)+P(i+1))/2;
    
    [r,alpha,mu_i,Cp] = Recuperacio (T_i, T_wall, P_i, V_i, R(i), Delta);
    f_i = Rugositat(V_i, R(i), epsilon, Rho_i, mu_i); %Skin friction coefficient
    
    Fregament = Fregament + f_i*Rho_i*V_i^2*A_i(i)/2;   
end

%Càlcul de la pressió normal al conducte


Momentum_out = m_in*V(Nvc+1)+P(Nvc+1)*S(Nvc+1);
Err_mom(n,m) = abs(Momentum_in-Fregament-Momentum_out)*100/Momentum_in;
if Err_mom(n,m)<0.1
    disp('Conservació momentum OKEY')
end

% 4.3 Conservació energia
%Energia entrada
[~,~,~,Cp_in] = Recuperacio(T_in,T_wall,P_in,V_in,R(1),Delta);
E_in = m_in*(Cp_in*T_in+(V_in^2)/2); 

%Energia sortida
[~,~,~,Cp_out] = Recuperacio(T(Nvc+1),T_wall,P(Nvc+1),V(Nvc+1),R(Nvc+1),Delta);
E_out = m_in*(Cp_out*T(Nvc+1)+(V(Nvc+1)^2)/2); 

%Càlcul del calor
Q = 0;
if Adiabatic~=1 
    for i=1:Nvc
        V_i = (V(i)+V(i+1))/2;
        T_i = (T(i)+T(i+1))/2;
        P_i = (P(i)+P(i+1))/2;

        [r, alpha, mu, Cp_i] = Recuperacio (T_i, T_wall, P_i, V_i, R(i), Delta); %Coeficients
        Tr = T_i + r*V_i^2/(Cp_i*2); %Temperatura de recuperació
        Q = Q + alpha*(T_wall-Tr)*Per(i)*Delta_z;   
    end
end

Err_ene(n,m) = abs(E_in-E_out+Q)*100/E_in;
if Err_ene(n,m)<0.1
    disp('Conservació energia OKEY')
end



% temps(n) = toc;%%% PER LESTUDI DE DELTA
% ERRORr(n) = abs(Rho_error-Rho(2001))*100/Rho_error ;%%% PER LESTUDI DE DELTA
% ERRORt(n) = abs(T_error-T(2001))*100/T_error;%%% PER LESTUDI DE DELTA
% ERRORv(n) = abs(V_error-V(2001))*100/V_error;%%% PER LESTUDI DE DELTA

% end

%Ara tenim els valors per un mach concret, ho guardem en una matriu que
%tingui tots els machs:
Vm(m,:) = V(1,:);
Pm(m,:) = P(1,:);
Tm(m,:) = T(1,:);
Rhom(m,:) = Rho(1,:);
Machm(m,:) = Mach(1,:);
sm(m,:) = s(1,:);
Sgenm(m,:) = Sgen(1,:);



end

%% 5. Impressió de resultats ---------------------------------------------


% 5.1. Impressio de resultats fisics
Impressio (z,Vm,Pm,Tm,Rhom,Machm,sm,Sgenm,M);

% 5.2. Impressió d'errors 
% Errors(N,Err_mom, Err_mass,Err_mom);

% 5.3. Impressió de convergència nodes
% Numeric( N,J_mean )








