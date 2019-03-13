function [alfa, r, f] = compressible(Ti, Tti, vi, Pi, Di)
dif2 = 30;
Trs = Ti; %Suposem Temp. recuperaci� = Temp. entrada
while dif2>10e-1
    
    Tref = Ti+0.22*(Trs-Ti); %Trs anir� canviant a cada iteraci�;
    rho = Pi/(287*Tref);
    
    %Propiertats del fluid al VC
    
    Cpi = 1022 - 0.1626*Tref+3.5025e-4*Tref^2;
    lambda = 3.807e-3+7.4e-5*Tref;
    mu = (2.53928e-5*sqrt(Tref/273.15))/(1+(122/Tref));
    
    %Grups adimensionals
    
    Re = rho*vi*Di/mu;
    Pr = mu*Cpi/lambda;
    
    %c�lcul del factor de relaxaci� i de Tr
    r = Pr^(1/3);
    Tr = Ti + r*vi^2/(2*Cpi); %Valor calculat de Tr.
    
    dif2 = abs(Tr-Trs); %�s |Trs-Tr|<error? 
    
    Trs = Tr; %Per la nova iteraci� suposem Trs = Tr. 
    
end
    %Ja s'ha calculat Tr amb �xit.
    
    Nu = 0.023*Re^0.8*Pr^0.4; %Sempre �s turbulent
    alfa = lambda*Nu/Di;
    f = 0.078 * Re^-0.2;
end
