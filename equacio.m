function [v1,v2,A_v,A_t,B_v,B_t,C_v,C_t] = equacio(Cpi,m_punt,alfa,r,f,rhoi,vi,Per,delta_x,S,P,T,v,Tti,i,j,w,R)

 %Resolem el sistema d'equacions

    A_v = m_punt + f*rhoi*abs(vi)*Per*delta_x/4;
    B_v = S;
    C_v = S*P(i,j,w)+(m_punt-f*rhoi*abs(vi)*Per*delta_x/4)*v(i,j,w);
    
    A_t = m_punt*Cpi + 0.5*alfa*Per*delta_x;
    B_t = 0.5*m_punt + (r*alfa*Per*delta_x)/(4*Cpi);
    C_t = (m_punt*Cpi - alfa*Per*delta_x*0.5)*T(i,j,w)+(0.5*m_punt-(r*alfa*Per*delta_x)/(4*Cpi))*v(i,j,w)^2+ alfa*Tti*Per*delta_x;
    
    A = A_v*A_t*S - B_v*B_t*m_punt*R;
    B = C_v*A_t*S;
    C = B_v*C_t*m_punt*R;

    
    %Ara tenim l'equació quadràtica
   
 
    det = B^2-4*A*C;
   
    if det < 0 %Discriminant negatiu, no te solució física
        error('El determinant és negatiu, revisar condicions inicials o geometria');
    end
    
    v1 = (B+sqrt(B^2-4*A*C))/(2*A); 
    v2 = (B-sqrt(B^2-4*A*C))/(2*A);
end

