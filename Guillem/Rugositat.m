function [ f ] = Rugositat(V,R,Er,rho,mu)
    Re = rho*V*R*2/mu; %Numero de Reynolds
    
    %Fem servir l'expressió de Churchill
    A = (2.457*log(1/((7/Re)^0.9 + 0.27*Er)))^16;
    B = (37530/Re)^16;
    f = 2*((8/Re)^12 + 1/(A+B)^(3/2))^(1/12);

  % f = 0.0625/((log10(Er/3.7 + 5.74/Re^0.9))^2);

   % f = 0.078*Re^-0.2;
end

