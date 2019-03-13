function [Tti] = conduccio(alfa,i,delta_x,lambdat,e,Tr,alfaext,Text,Per,variable,j,w,Treferencia,Tt,N)

    
    if variable == true
    ae = (lambdat*e^2*pi)/delta_x; %area aproximada de pi*radi^2*espessor.
    aw = (lambdat*e^2*pi)/delta_x;
    bp = alfa*Tr*Per*delta_x+alfaext*Text*Per*delta_x;
    ap = ae+aw+alfa*Per*delta_x+alfaext*Per*delta_x; %considerem superficie interior i exterior aprox. igual.
    if i == 1
        ap = 1;
        aw = 0;
        bp= Treferencia;
        ae = 0;
    end
    
    if i == 1
        Tti = (ae*Tt(i+1,j,w)+bp)/ap; %nomes conduccio per la dreta
        elseif i == N   
        ap = aw+alfa*Per*delta_x+alfaext*Per*delta_x;
        Tti = (aw*Tt(i-1,j,w)+bp)/ap; %Nomes conduccio per l'esquerra.
         else
        Tti = (ae*Tt(i+1,j,w)+aw*Tt(i-1,j,w)+bp)/ap;
    end

    
    else
        Tti = 300;
    end
end
