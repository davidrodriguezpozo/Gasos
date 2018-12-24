function [ r, alpha ] = Recuperacio (T_inf, T_wall, P_inf, v, R, Delta)
    Tr = T_inf; %Suposem la temperatura de recuperació
    for i=1:100 %Fixem un màxim de 100 iteracions
        T_ref = (T_wall+T_inf)/2 + 0.22*(Tr-T_inf); %Temperatura de referencia
        
        %Calculem propietats termofísique amb aquesta T_ref:
        rho = P_inf / (287*T_ref); %Densitat
        if T_ref>=1500
            mu = 2.5393*10^(-5) * sqrt(T_ref/273.15) / (1+(122/T_ref)); %Viscositat dinàmica
        end
        if T_ref<1500
            mu = (1.458*10^(-6)*T_ref^1.5)/(T_ref+110.4);
        end            
        lambda = 2.648*10^(-3)*sqrt(T_ref)/(1+(245.4/T_ref)*10^(-12/T); %Conductivitat tèrmica
        Cp = 1034.09-2.849*10^(-1)*T_ref+7.817*10^(-4)*T_ref^2-4.971*10^(-7)*T_ref^3+1.088*10^(-10)*T_ref^4;
        
        %Re i Pr local
        Re = rho*v*2*R/mu;
        Pr = mu*Cp/lambda;
        
        %Laminar o turbulent?
        if Re<2000 %Laminar
            r = Pr^(1/2);
        end
        if Re>2000 %Turbulent
            r = Pr^(1/3);
        end
        
        Tr_calc = T_inf + r*v^2/(Cp*2); %Calculem la nova temperatura de recuperació
        
        %Evaluem convergencia
        if abs(Tr-Tr_calc)<Delta
            break
        end
        Tr = Tr_calc;
    end
    
    %Laminar o turbulent?
    if Re<2000 %Laminar
       Nu = 3.66;
    end
    
    if Re>2000 %Turbulent
       Nu = 0.023*Re^0.8*Pr^0.4;
    end
    
    alpha = Nu*lambda/(2*R);    
end

