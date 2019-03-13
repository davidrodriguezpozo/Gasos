%----------DINÀMICA DE GASOS-------------

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%RESOLUCIÓ ANALÍTICA DE PROBLEMA COMBINAT. 

%TUB AMB FLUID INTERIOR I CONVECCIÓ EXTERIOR.

clear all;

%%-----------CÀLCULS PREVIS---------------

%%Entrada de dades per triar el tipus problema Supersonic-Subsonic-Ones de
%%xoc

onesdexoc = false; %posar true si volem ones de xoc.
Quantonesdexoc = 1; %sleccionar quantes ones de xoc es volen, si es volen ones de xoc, posar 3 (predeterminat)
M = [0.3;0.6;0.7;0.78];  %Mach d'entrada a M > 0.78 es torna sònic. 
adiabatic = false; %En cas de que sigui adiabàtic --> true
variable = true; %Posem false si la temperatura del tub és constant. 

%%Dades geomètriques del tub 

Di = 0.01; %Diàmetre del tub [m]
ri = Di/2; %Radi [m]
L = 0.05; %Longitud del tub [m]
epsilon = 0.004; 
S = pi*ri^2; %Superfície del VC [m^2]
Per = Di*pi; %Perímetre del VC. [m^2]
e = 0.001; %Espessor del tub d'1mm [m]


%%Dades del fluid i del tub

lambdat = 300; %Coeficient de conductivitat termica del tub (posem 300)
Treferencia = 300; %Posem una temperatura per iniciar el tub. [K]
R = 287; %Constant dels gasos per l'aire
gamma = 1.4; %Constant adiab. de l'aire

Tin = 400; %Temperatura d'entrada [K]
pin = 5e5; %Pressió d'entrada [Pa]
rhoin = pin/(287*Tin); %Densitat a l'entrada [Kg/m^3]


onadexoc = false; %Posem false al principi per evitar una ona de xoc.
alfaext = 100; %alfa exterior, considerem 100. 
Text = 20+273; %Temperatura de l'aire exterior (20ºC)
m = Quantonesdexoc;

%Discretització del tub

N = 1000; %numero de VC
delta_x = L/N; %Longitud del VC
x = 0:delta_x:L; %Vector de la longitud del tub 

%Definicio de tots els vectors (NxMxm) (VC,Mach d'entrada,Ona de xoc)

Tt = zeros(N,length(M),m);
Tt_s = zeros(N,length(M),m);
T_s = zeros(N+1,length(M),m);
V_s = zeros(N+1,length(M),m);
P_s = zeros(N+1,length(M),m);
rho_s = zeros(N+1,length(M),m);
v = zeros(N+1,length(M),m);
T = zeros(N+1,length(M),m);
P = zeros(N+1,length(M),m);
rho = zeros(N+1,length(M),m);
q = zeros(N,length(M),m);
T_r = zeros(N,length(M),m);
Sgen = zeros(N+1,length(M),m);
Sgen1_v = zeros(N,length(M),m);
Sgen2_v = zeros(N,length(M),m);
S_spec = zeros(N+1,length(M),m);
det_v = zeros(N,length(M),m);
P1_v = zeros(N,length(M),m);
P2_v = zeros(N,length(M),m);
T1_v = zeros(N,length(M),m);
T2_v = zeros(N,length(M),m);
v1_v = zeros(N,length(M),m);
v2_v = zeros(N,length(M),m);
Mach = zeros(N+1,length(M),m);
f_v = zeros(N,length(M),m);
m_v = zeros(length(M),m);
vin_v = zeros(length(M),m);



%% Resolució del problema
for w = 1:m %Per guardar cada resolució en diferents vector per ones de xoc. 

for j = 1:length(M) %Variem el Mach a l'entrada
    Min = M(j);
    vin = sqrt(gamma*R*Tin)*Min; %Velocitat d'entrada 
    m_punt = S*vin*rhoin; %Cabal màssic [Kg/s]
    vin_v(j,w) = vin;
    m_v(j,w) = m_punt;
    
    %%Omplim els vectos suposats amb cariables conegudes, en aquest cas les
    %%d'entrada
    
    for i = 1:N
    Tt(i,j,w) = Treferencia;
    Tt_s(i,j,w) = Treferencia;    
    end
    
    for i = 1:N+1
    T_s(i,j,w) = Tin;
    V_s(i,j,w) = vin;
    P_s(i,j,w) = pin;
    rho_s(i,j,w) = rhoin;
    T(i,j,w) = Tin;
    v(i,j,w) = vin;
    P(i,j,w) = pin;
    rho(i,j,w) = rhoin;
    Mach(i,j,w) = Min;
    end
    
for i = 1:N
    
    T_s(i+1,j,w) = T(i,j,w); %Suposem variables de sortida = entrada.
    P_s(i+1,j,w) = P(i,j,w);
    V_s(i+1,j,w) = v(i,j,w);
    rho_s(i+1,j,w) = rho(i,j,w);
    if i>1
    Tt_s(i,j,w) = Tt(i-1,j,w);
    end
    
    dif = 10e10; %Valor arbitrari per poder iterar un altre cop
    y=0;
    while dif > 1e-10 && y<1000 %Evitem que faci més de 1000 iteracions.

    y=y+1;
    Ti = 0.5*(T(i,j,w)+T(i+1,j,w)); % Variables mitges del VC.
    vi = 0.5*(v(i,j,w)+v(i+1,j,w)); 
    Pi = 0.5*(P(i,j,w)+P(i+1,j,w));
    rhoi = Pi/(287*Ti);
    Tti = Tt(i,j,w);
    %Càlcul dels coeficient alfa, r i f
    
    [alfa, r, f] = compressible(Ti, Tti, vi, Pi, Di); 
    
    %Un cop calculats seguim amb el càlcul de Tr
    f_v(i,j,w) = f;
    Cpi = 1034.09-2.849*10^(-1)*Ti+7.817*10^(-4)*Ti^2-4.971*10^(-7)*Ti^3+1.088*10^(-10)*Ti^4;
    Tr = Ti + r*vi^2/(2*Cpi);
    
    %Calcul de conduccio
    
    [Tti] = conduccio(alfa,i,delta_x,lambdat,e,Tr,alfaext,Text,Per,variable,j,w,Treferencia,Tt,N);
    
    Tt(i,j,w) = Tti;
    
    q(i,j,w) = alfa*(Tti-Tr)*Per*delta_x; %podem calcular el flux de calor amb Tr;
    
    [v1,v2,A_v,A_t,B_v,B_t,C_v,C_t] = equacio(Cpi,m_punt,alfa,r,f,rhoi,vi,Per,delta_x,S,P,T,v,Tti,i,j,w,R);
   
    v1_v(i,j,w) = v1;
    v2_v(i,j,w) = v2;
    
    % Calculem totes les propietats amb les dues velocitats.
    
    Vol = S*delta_x;
    P1 = (C_v-A_v*v1)/B_v;
    T1 = (C_t-B_t*v1^2)/A_t;
    P2 = (C_v-A_v*v2)/B_v;
    T2 = (C_t-B_t*v2^2)/A_t;
    
    P1_v(i,j,w) = P1; %Guardem en vectors per poder veure-ho desrpés si cal
    P2_v(i,j,w) = P2;
    T2_v(i,j,w) = T2;
    T1_v(i,j,w) = T1;
    
    %Trobem l'entropia generada amb les dues velocitats.
    if adiabatic == true
        alfa = 0;
    end
    
    if onesdexoc == true
    if w == 1
        
        if x(i) == L/2 
        onadexoc = true;
        end
        
    elseif w == 2
        if x(i) == 4*L/5
            onadexoc = true;
        end
        
    else
           if x(i) == 3*L/4
                onadexoc = true;
           end
    end
    end
    
    
    
    Sgen1 = 1/Vol*(m_punt*(Cpi*log(T1/T(i,j,w))-R*log(P1/P(i,j,w))) - alfa*(Tti-Tr)*Per*delta_x/Tti);
    Sgen2 = 1/Vol*(m_punt*(Cpi*log(T2/T(i,j,w))-R*log(P2/P(i,j,w))) - alfa*(Tti-Tr)*Per*delta_x/Tti);
    
    Sgen1_v(i,j,w) = Sgen1; %tornem a guardar per si acàs. 
    Sgen2_v(i,j,w) = Sgen2;
    
    
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
        v(i+1,j,w) = v2;
        T(i+1,j,w) = T2;
        P(i+1,j,w) = P2;
        Sgen(i+1,j,w) = Sgen2;
        %disp('Entropia 1 és imaginaria, entropia 2 positiva');
    end
    if isreal(Sgen2) == 0 && Sgen1 > 0 %Comprovem si Sgen2 es real o no
        v(i+1,j,w) = v1;
        T(i+1,j,w) = T1;
        P(i+1,j,w) = P1;
        Sgen(i+1,j,w) = Sgen1;
        disp('Entropia 2 és imaginaria, entropia 1 positiva');
    end
    
    
    if abs(Sgen1) > abs(Sgen2) %Cas en que les dues són positives, ens quedem amb la segona velocitat
        v(i+1,j,w) = v2;
        T(i+1,j,w) = T2;
        P(i+1,j,w) = P2;
        Sgen(i+1,j,w) = Sgen2;
    else
        v(i+1,j,w) = v1;
        T(i+1,j,w) = T1;
        P(i+1,j,w) = P1;
        Sgen(i+1,j,w) = Sgen1;
    end
    
    if onadexoc == true %Si hi ha ona de xoc ens quedem amb la que genera més entropia. 
         if abs(Sgen1) > abs(Sgen2) %Cas en que les dues són positives, ens quedem amb la segona velocitat
        v(i+1,j,w) = v1;
        T(i+1,j,w) = T1;
        P(i+1,j,w) = P1;
        Sgen(i+1,j,w) = Sgen1;
    else
        v(i+1,j,w) = v2;
        T(i+1,j,w) = T2;
        P(i+1,j,w) = P2;
        Sgen(i+1,j,w) = Sgen2;
         end
         onadexoc = false; %Perque no torni a entrar
    end
    
    
    rho(i+1,j,w) = P(i+1,j,w)/(R*T(i+1,j,w));
    c = sqrt(gamma*R*T(i,j,w));
    Mach(i+1,j,w) = v(i+1,j,w)/c;
    
    %Guardem i comprovem les diferències entre suposat-calculat.
    
    difvector(1) = abs(v(i+1,j,w)-V_s(i+1,j,w));
    difvector(2) = abs(T(i+1,j,w)-T_s(i+1,j,w));
    difvector(3) = abs(P(i+1,j,w)-P_s(i+1,j,w));
    difvector(4) = abs(Tt(i,j,w)-Tt_s(i,j,w));
    
    dif = max(difvector);

    %Canviem el vector suposat pel calculat i tornem a iterar.
    
    V_s(i+1,j,w) = v(i+1,j,w);
    P_s(i+1,j,w) = P(i+1,j,w);
    T_s(i+1,j,w) = T(i+1,j,w);  
    Tt_s(i,j,w) = Tt(i,j,w);
    if i ==1
        S_spec(i,j,w) = 0;
    else
        S_spec(i,j,w) = 0;
    S_spec(i,j,w) = S_spec(i-1,j,w)+Sgen(i,j,w)*delta_x*S*(1/m_punt)+q(i,j,w)/Tt(i,j,w);
    end
    
    end
    
end

end
end



%% Verificació del codi

for m = 1:w
for j =1:length(M)
vin = vin_v(j,w);
m_punt = m_v(j,w);  

% Massa

massa_out = v(N+1,j,w)*rho(N+1,j,w)*S;
dif_mas = abs(m_v(j,w)-massa_out)*100/m_v(j,w);
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
    freg = freg + f_v(i,j,w)*rhoi*vi^2*pi*Di*delta_x/2; %Sumem la força de fregament total
end

Momentumout = m_punt*v(N+1,j,w)+P(N+1,j,w)*S;
dif_mom = abs(Momentumin-freg-Momentumout)*100/Momentumin;
if dif_mom<0.1
    disp('Conservacio momentum correcta per Mach =')
    disp('')
    disp(M(j))
    disp('')
end


% Energia

%Energia entrada
Cpin = 1022-0.1626*Tin+3.5025*10^(-4)*Tin^2;
Ein = m_punt*(Cpin*Tin+(vin^2)/2); 

%Energia sortida
Cpout = 1022-0.1626*T(N+1,j,w)+3.5025*10^(-4)*T(N+1,j,w)^2;
Eout = m_punt*(Cpout*T(N+1,j,w)+(v(N+1,j,w)^2)/2); 

%E = Ein - Q; Hem de calcular Q (calor).
 Q = 0;
for i=1:N
        vi = (v(i,j,w)+v(i+1,j,w))/2;
        Ti = (T(i,j,w)+T(i+1,j,w))/2;
        Pi = (P(i,j,w)+P(i+1,j,w))/2;
        Tti = Tt(i,j,w);
        Cpi = 1022 - 0.1626*Ti+3.5025e-4*Ti^2;
        [alfa, r, f] = compressible(Ti, Tti, vi, Pi, Di);
        Tr = Ti + r*vi^2/(Cpi*2); %Temperatura de recuperació
        Q = Q + alfa*(Tti-Tr)*pi*Di*delta_x; 
end


dif_ene = abs(Ein-Eout+Q)*100/Ein;
if dif_ene<0.1
    disp('Conservacio energia correcta per Mach = ') 
    disp('')
    disp(M(j))
    disp(' ')
end


end
end





%% Plots per les diferents variables

grafiques(P,T,rho,Mach,v,Sgen,S_spec,Tt,x,m,M,N)

