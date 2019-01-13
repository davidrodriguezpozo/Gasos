%----------DINÀMICA DE GASOS-------------

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%RESOLUCIÓ ANALÍTICA DE PROBLEMA COMBINAT. 

%TUB AMB FLUID INTERIOR I CONVECCIÓ EXTERIOR.

clear all;

%%-----------CÀLCULS PREVIS---------------



%%Dades geomètriques

Di = 0.01; %Diàmetre del tub [m]
ri = Di/2; %[m]
L = 0.05; %Longitud del tub [m]
epsilon = 0.004; 
S = pi*ri^2; %Superfície del VC [m^2]
Per = Di*pi; %Perímetre del VC. [m^2]

%%Dades del fluid

R = 287; %Constant dels gasos per l'aire
gamma = 1.4; %Constant adiab. de l'aire
M = [0.3;0.6;0.7;0.7851];  %Mach d'entrada a M > 0.7851 es torna sònic. 
Tin = 400; %Temperatura d'entrada [K]
pin = 5e5; %Pressió d'entrada [Pa]
rhoin = pin/(287*Tin); %Densitat a l'entrada [Kg/m^3]
Tt = 300; %Temperatura de la paret del tub [K]
adiabatic = false; %En cas de que sigui adiabàtic --> true


%Discretització del tub

N = 1000; %numero de VC
delta_x = L/N; %Longitud del VC
x = 0:delta_x:L; %Vector de la longitud del tub 

%Definicio de tots els vectors (NxM)

Tt = zeros(N,length(M));
ap = zeros(N,length(M));
aw = zeros(N,length(M));
ae = zeros(N,length(M));
bp = zeros(N,length(M));
T_s = zeros(N+1,length(M));
V_s = zeros(N+1,length(M));
P_s = zeros(N+1,length(M));
rho_s = zeros(N+1,length(M));
v = zeros(N+1,length(M));
T = zeros(N+1,length(M));
P = zeros(N+1,length(M));
rho = zeros(N+1,length(M));
q = zeros(N,length(M));
T_r = zeros(N,length(M));
Sgen = zeros(N+1,length(M));
Sgen1_v = zeros(N,length(M));
Sgen2_v = zeros(N,length(M));
det_v = zeros(N,length(M));
P1_v = zeros(N,length(M));
P2_v = zeros(N,length(M));
T1_v = zeros(N,length(M));
T2_v = zeros(N,length(M));
v1_v = zeros(N,length(M));
v2_v = zeros(N,length(M));
Mach = zeros(N+1,length(M));
f_v = zeros(N,length(M));
m_v = zeros(length(M));
vin_v = zeros(length(M));
Trs = 300; %Suposem la temperatura de recuperació inicial
%Omplim els vectors amb les variables suposades.




%% Resolució del problema

for j = 1:length(M) %Variem el Mach a l'entrada
    Min = M(j);
    vin = sqrt(gamma*R*Tin)*Min; %Velocitat d'entrada 
    m_punt = S*vin*rhoin; %Cabal màssic [Kg/s]
    vin_v(j) = vin;
    m_v(j) = m_punt;
    
    for i = 1:N+1
    T_s(i,j) = Tin;
    V_s(i,j) = vin;
    P_s(i,j) = pin;
    rho_s(i,j) = rhoin;
    T(i,j) = Tin;
    v(i,j) = vin;
    P(i,j) = pin;
    rho(i,j) = rhoin;
    Mach(i,j) = Min;
    end
    
for i = 1:N

    T_s(i+1,j) = T(i,j); %Suposem variables de sortida = entrada.
    P_s(i+1,j) = P(i,j);
    V_s(i+1,j) = v(i,j);
    rho_s(i+1,j) = rho(i,j);
    dif = 10e10; %Valor arbitrari per poder iterar un altre cop
    y=0;
    while dif > 1e-10 && y<1000 %Evitem que faci més de 1000 iteracions.

    y = y+1; %Comptador d'iteracions
    
    Ti = 0.5*(T(i,j)+T(i+1,j)); % Variables mitges del VC.
    vi = 0.5*(v(i,j)+v(i+1,j)); 
    Pi = 0.5*(P(i,j)+P(i+1,j));
    rhoi = Pi/(287*Ti);
    
    %Càlcul dels coeficient alfa, r i f
    
    [alfa, r, f] = compressible(Ti, Tt, vi, Pi, Di); 
    
    %Un cop calculats seguim amb el càlcul de Tr
    f_v(i,j) = f;
    Cpi = 1034.09-2.849*10^(-1)*Ti+7.817*10^(-4)*Ti^2-4.971*10^(-7)*Ti^3+1.088*10^(-10)*Ti^4;
    Tr = T(i) + r*vi^2/(2*Cpi);
    q(i,j) = alfa*(Tt-Tr)*Per*delta_x; %podem calcular el flux de calor amb Tr;
    
    %Resolem el sistema d'equacions

    A_v = m_punt + f*rhoi*abs(vi)*Per*delta_x/4;
    B_v = S;
    C_v = S*P(i,j)+(m_punt-f*rhoi*abs(vi)*Per*delta_x/4)*v(i,j);
    
    A_t = m_punt*Cpi + 0.5*alfa*Per*delta_x;
    B_t = 0.5*m_punt + (r*alfa*Per*delta_x)/(4*Cpi);
    C_t = (m_punt*Cpi - alfa*Per*delta_x*0.5)*T(i,j)+(0.5*m_punt-(r*alfa*Per*delta_x)/(4*Cpi))*v(i,j)^2+ alfa*Tt*Per*delta_x;
    
    A = A_v*A_t*S - B_v*B_t*m_punt*R;
    B = C_v*A_t*S;
    C = B_v*C_t*m_punt*R;

    
    %Ara tenim l'equació quadràtica
   
    
    det = B^2-4*A*C;
    det_v(i,j) = det;
    
    if det < 0 %Discriminant negatiu, no te solució física
        error('El determinant és negatiu');
    end
    
    v1 = (B+sqrt(B^2-4*A*C))/(2*A); 
    v2 = (B-sqrt(B^2-4*A*C))/(2*A);
    v1_v(i,j) = v1;
    v2_v(i,j) = v2;
    
    % Calculem totes les propietats amb les dues velocitats.
    
    Vol = S*delta_x;
    P1 = (C_v-A_v*v1)/B_v;
    T1 = (C_t-B_t*v1^2)/A_t;
    P2 = (C_v-A_v*v2)/B_v;
    T2 = (C_t-B_t*v2^2)/A_t;
    
    P1_v(i,j) = P1; %Guardem en vectors per poder veure-ho desrpés si cal
    P2_v(i,j) = P2;
    T2_v(i,j) = T2;
    T1_v(i,j) = T1;
    
    %Trobem l'entropia generada amb les dues velocitats.
    if adiabatic == true
        alfa = 0;
    end
    
    Sgen1 = 1/Vol*(m_punt*(Cpi*log(T1/T(i,j))-R*log(P1/P(i,j))) - alfa*(Tt-Tr)*Per*delta_x/Tt);
    Sgen2 = 1/Vol*(m_punt*(Cpi*log(T2/T(i,j))-R*log(P2/P(i,j))) - alfa*(Tt-Tr)*Per*delta_x/Tt);
    
    Sgen1_v(i,j) = Sgen1; %tornem a guardar per si acàs. 
    Sgen2_v(i,j) = Sgen2;
    
    if Sgen1 < 0 && Sgen2 < 0  %Si les dues són negatives, no té solució.
        error('Les dues entropies són negatives');
    end
    
    if isreal(Sgen1) == 0 && isreal(Sgen2) == 0 %Fem que si hi ha una netropia im. digui quina és.
        error('Les dues entropies són imaginaries');
    end
    
    if (isreal(Sgen1) == 0 && Sgen2 < 0) || (isreal(Sgen2) == 0 && Sgen1 < 0)
        disp('Entropies no es poden resoldre');
    end
        
    if isreal(Sgen1) == 0 && Sgen2 > 0 %Comprovem si Sgen1 es real o no
        v(i+1,j) = v2;
        T(i+1,j) = T2;
        P(i+1,j) = P2;
        Sgen(i+1,j) = Sgen2;
        %disp('Entropia 1 és imaginaria, entropia 2 positiva');
    end
    if isreal(Sgen2) == 0 && Sgen1 > 0 %Comprovem si Sgen2 es real o no
        v(i+1,j) = v1;
        T(i+1,j) = T1;
        P(i+1,j) = P1;
        Sgen(i+1,j) = Sgen1;
        disp('Entropia 2 és imaginaria, entropia 1 positiva');
    end
    if abs(Sgen1) > abs(Sgen2) %Cas en que les dues són positives
        v(i+1,j) = v2;
        T(i+1,j) = T2;
        P(i+1,j) = P2;
        Sgen(i+1,j) = Sgen2;
    else
        v(i+1,j) = v1;
        T(i+1,j) = T1;
        P(i+1,j) = P1;
        Sgen(i+1,j) = Sgen1;
    end
    
    rho(i+1,j) = P(i+1,j)/(R*T(i+1,j));
    c = sqrt(gamma*R*T(i,j));
    Mach(i+1,j) = v(i+1,j)/c;
    
    %Guardem i comprovem les diferències entre suposat-calculat.
    
    difvector(1) = abs(v(i+1,j)-V_s(i+1,j));
    difvector(2) = abs(T(i+1,j)-T_s(i+1,j));
    difvector(3) = abs(P(i+1,j)-P_s(i+1,j));
    
    dif = max(difvector);

    %Canviem el vector suposat pel calculat i tornem a iterar.
    
    V_s(i+1,j) = v(i+1,j);
    P_s(i+1,j) = P(i+1,j);
    T_s(i+1,j) = T(i+1,j);      
    
   end
    
end

end




%% Verificació del codi


for j =1:length(M)
vin = vin_v(j);
m_punt = m_v(j);  

% Massa

massa_out = v(N+1,j)*rho(N+1,j)*S;
dif_mas = abs(m_v(j)-massa_out)*100/m_v(j);
if dif_mas<0.1
    disp('Conservacio massa correcta per Mach = ')
    disp('')
    disp(M(j))
    disp('')
end

% Momentum

Momentumin = m_punt*vin+pin*S; %Calculem el momentum a l'entrada
freg = 0;
for i=1:N
    vi = (v(i,j)+v(i+1,j))/2;
    rhoi = (rho(i,j)+rho(i+1,j))/2;
    freg = freg + f_v(i,j)*rhoi*vi^2*pi*Di*delta_x/2; %Sumem la força de fregament total
end

Momentumout = m_punt*v(N+1,j)+P(N+1,j)*S;
dif_mom = abs(Momentumin-freg-Momentumout)*100/Momentumin;
if dif_mom<0.1
    disp('Conservacio momentum correcta per Mach =')
    disp('')
    disp(M(j))
    disp('')
end


% 4.3 Conservació energia

%Energia entrada
Cpin = 1022-0.1626*Tin+3.5025*10^(-4)*Tin^2;
Ein = m_punt*(Cpin*Tin+(vin^2)/2); 

%Energia sortida
Cpout = 1022-0.1626*T(N+1,j)+3.5025*10^(-4)*T(N+1,j)^2;
Eout = m_punt*(Cpout*T(N+1,j)+(v(N+1,j)^2)/2); 

%Calculem les pèrdues per calor
 Q = 0;
for i=1:N
        vi = (v(i,j)+v(i+1,j))/2;
        Ti = (T(i,j)+T(i+1,j))/2;
        Pi = (P(i,j)+P(i+1,j))/2;
        Cpi = 1034.09-2.849*10^(-1)*Ti+7.817*10^(-4)*Ti^2-4.971*10^(-7)*Ti^3+1.088*10^(-10)*Ti^4;
        [alfa, r, f] = compressible(Ti, Tt, vi, Pi, Di);
        Tr = Ti + r*vi^2/(Cpi*2); %Temperatura de recuperació
        Q = Q + alfa*(Tt-Tr)*pi*Di*delta_x;   
end


dif_ene = abs(Ein-Eout+Q)*100/Ein;
if dif_ene<0.1
    disp('Conservacio energia correcta per Mach = ') 
    disp('')
    disp(M(j))
    disp(' ')
end


end


%% Plots per les diferents variables
figure;

subplot(2,2,1);
plot(x,P);
grid on;
title('Pressi\''o al llarg del tub');
xlabel('Longitud del tub [m]');
ylabel('Pressi\''o [Pa]');
ylim([300000 550000]);


subplot(2,2,2);
plot(x,T);
grid on;
title('Temperatura');
xlabel('Longitud del tub [m]');
ylabel('Temperatura [K]')
ylim([370 405]);

subplot(2,2,3);
plot(x,v);
grid on;
title('Velocitat');
xlabel('Longitud del tub [m]');
ylabel('Velocitat [$\frac{m}{s}$]')
ylim([100 400]);

subplot(2,2,4);
plot(x,Mach);
grid on;
title('Mach');
legendCell = cellstr(num2str(M, 'M=%-d'));
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Nombre de Mach');
ylim([0.2 1]);


figure;

subplot(2,1,1);
plot(x,rho);
grid on;
title('Densitat al llarg del tub');
xlabel('Longitud del tub [m]');
ylabel('Densitat [$\frac{Kg}{m^3}$]');
ylim([3.5 4.5]);

subplot(2,1,2);
plot(x,Sgen);
grid on;
title('Entropia generada');
legendCell = cellstr(num2str(M, 'M=%-d'));
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Entropia generada \textit{Sgen} [$\frac{J}{K}$]');

