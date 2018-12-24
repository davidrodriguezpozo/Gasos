function [alfa, r, f, Cpi] = compressible(Ti, Tt, vi, Pi, Di)
dif2 = 30;
Trs = Ti; %Suposem Temp. recuperaci� = Temp. entrada
while dif2>10e-10
    
    Tref = 0.5*(Tt+Ti)+0.22*(Trs-Ti); %Trs anir� canviant a cada iteraci�;
    rho = Pi/(287*Tref);
    
    %Propiertats del fluid al VC
    
    mu = (2.5393e-5*sqrt(Tref/273.15))/(1+122/Tref);
    lambda = 3.807e-3+7.4e-5*Tref;
    Cpi = 1022 - 0.166*Tref + 3.5025e-4*Tref^2;
    
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
    C = 0.023; m=0.8; n=0.4; K = 1;
    Nu = C*Re^m*Pr^n*K;
    alfa = lambda*Nu/Di;
    f = 0.078 * Re^-0.2;
end
