function [alfa, r, f] = compressible(Ti, Tt, vi, Pi, Di)
dif2 = 30;
Trs = Ti; %Suposem Temp. recuperació = Temp. entrada
while dif2>10e-1
    
    Tref = 0.5*(Tt+Ti)+0.22*(Trs-Ti); %Trs anirà canviant a cada iteració;
    rho = Pi/(287*Tref);
    
    %Propiertats del fluid al VC
    
        %mu = (2.5393e-5*sqrt(Tref/273.15))/(1+122/Tref);
        %lambda = 3.807e-3+7.4e-5*Tref;
        %Cpi = 1022 - 0.166*Tref + 3.5025e-4*Tref^2;
    
    mu = (1.458*10^(-6)*Tref^1.5)/(Tref+110.4);
    lambda = 2.648*10^(-3)*sqrt(Tref)/(1+(245.4/Tref)*10^(-12/Tref));
    Cpi = 1034.09-2.849*10^(-1)*Tref+7.817*10^(-4)*Tref^2-4.971*10^(-7)*Tref^3+1.088*10^(-10)*Tref^4;
    
    
    %Grups adimensionals
    
    Re = rho*vi*Di/mu;
    Pr = mu*Cpi/lambda;
    
    %càlcul del factor de relaxació i de Tr
    r = Pr^(1/3);
    Tr = Ti + r*vi^2/(2*Cpi); %Valor calculat de Tr.
    
    dif2 = abs(Tr-Trs); %És |Trs-Tr|<error? 
    
    Trs = Tr; %Per la nova iteració suposem Trs = Tr. 
    
end
    %Ja s'ha calculat Tr amb èxit.
    
    Nu = 0.023*Re^0.8*Pr^0.4; %Sempre és turbulent
    alfa = lambda*Nu/Di;
    f = 0.078 * Re^-0.2;
end
