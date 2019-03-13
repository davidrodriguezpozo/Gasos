function grafiques(P,T,rho,Mach,v,Sgen,S_spec,Tt,x,m,M,N)


figure;


subplot(1,2,1);
for i = 1:m
plot(x,P(:,:,i)); hold on;
end
grid on;
title('Pressi\''o');
xlabel('Longitud del tub [m]');
ylabel('Pressi\''o [Pa]');


subplot(1,2,2);
for i = 1:m
plot(x,T(:,:,i)); hold on;
end
grid on;
title('Temperatura');
legendCell = cellstr(num2str(M, 'M=%-d'));
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Temperatura [K]')



figure;

subplot(1,2,1);
for i=1:m
plot(x,v(:,:,i)); hold on;
end
grid on;
title('Velocitat');
xlabel('Longitud del tub [m]');
ylabel('Velocitat [$\frac{m}{s}$]')

subplot(1,2,2);
for i=1:m
plot(x,Mach(:,:,i)); hold on;
end
grid on;
title('Mach');
legendCell = cellstr(num2str(M, 'M=%-d'));
format shortG;
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Nombre de Mach');





figure;

subplot(1,2,1);
S_spec(N+1,:,:) = S_spec(N,:,:);
for i=1:m
plot(x,S_spec(:,:,i)); hold on;
end
grid on;
title('Entropia especifica');
xlabel('Longitud del tub [m]');
ylabel('Entropia [$\frac{Kg}{m^3}$]');

subplot(1,2,2);
for i=1:m
plot(x,Sgen(:,:,i)); hold on;
end
grid on;
title('Entropia generada');
legendCell = cellstr(num2str(M, 'M=%-d'));
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Entropia generada \textit{Sgen} [$\frac{J}{K}$]');

figure;
plot(Tt(:,:,1));
legendCell = cellstr(num2str(M, 'M=%-d'));
legend(legendCell);
title('Temperatura del tub $T_{tub}$');
xlabel('Longitud del tub [m]');
ylabel('Temperatura [K]');


figure;
for i=1:m
plot(x,rho(:,:,i)); hold on;
end
grid on;
title('Densitat');
legendCell = cellstr(num2str(M, 'M=%-d'));
format shortG;
legend(legendCell);
xlabel('Longitud del tub [m]');
ylabel('Densitat [$\frac{Kg}{m^3}$]');

end
