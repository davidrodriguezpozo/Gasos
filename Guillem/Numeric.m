function [] = Numeric( N,J_mean,ERRORt,ERRORr,ERRORv,D,temps )


%NVC
plot(N,J_mean)
title('Iteracions per node')
xlabel('N�mero volums de control')
ylabel('Mitjana d''iteracions')
grid on

%Delta convergencia
ERR = ERRORt+ERRORr+ERRORv;
ERR = ERR/3;

figure()
title('Influ�ncia del factor de converg�ncia')
set(gca, 'XDir','reverse')
yyaxis left
loglog(D,ERR)
ylabel('Error [%]')
hold on
yyaxis right
plot(D,temps)
ylabel('Temps [s]')
xlabel('Factor de converg�ncia delta')
grid on
xlim([10^-20 0.1])

end

