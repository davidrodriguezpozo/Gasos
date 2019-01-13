function [] = Numeric( N,J_mean,ERRORt,ERRORr,ERRORv,D,temps )


%NVC
plot(N,J_mean)
title('Iteracions per node')
xlabel('Número volums de control')
ylabel('Mitjana d''iteracions')
grid on

%Delta convergencia
ERR = ERRORt+ERRORr+ERRORv;
ERR = ERR/3;

figure()
title('Influència del factor de convergència')
set(gca, 'XDir','reverse')
yyaxis left
loglog(D,ERR)
ylabel('Error [%]')
hold on
yyaxis right
plot(D,temps)
ylabel('Temps [s]')
xlabel('Factor de convergència delta')
grid on
xlim([10^-20 0.1])

end

