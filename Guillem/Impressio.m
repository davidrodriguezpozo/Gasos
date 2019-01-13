function [] = Impressio( z,V,P,T,Rho,Mach,s,Sgen,M)


figure()
subplot(1,2,1)
plot(z,V)
xlabel ('x [m]')
ylabel('v [m/s]')
title('Velocitat')
grid on
xlim([0,z(end)])
%ylim([40, inf])
hold on
for i=1:(length(M))
    legendInfo{i}=['M_{in}: ' num2str(M(i))];
    hold on
end
legend(legendInfo);


subplot(1,2,2)
plot(z,Mach)
title('Mach')
xlabel ('x [m]')
ylabel('Ma')
grid on
xlim([0,z(end)])
%ylim([0.1, 0.9])



% PRESSIÓ I TEMPERATURA --------------------------------------------------
figure()
subplot(1,2,1)
plot(z,P)
title('Pressió')
xlabel ('x [m]')
ylabel('P [Pa]')
grid on
legend(legendInfo);
xlim([0,z(end)])

subplot(1,2,2)
plot(z,T)
title('Temperatura')
xlabel ('x [m]')
ylabel('T [K]')
grid on
xlim([0,z(end)])








% DENSITAT I ENTROPIA -----------------------------------------
figure()
plot(z,Rho)
title('Densitat')
xlabel ('x [m]')
ylabel('Rho [kg/m^3]')
grid on
xlim([0,z(end)])
legend(legendInfo);




figure()
subplot(1,2,1)
plot(z,s)
title('Entropia acumulada')
xlabel ('x [m]')
ylabel('S [J/kgK]')
grid on
xlim([0,z(end)])
legend(legendInfo);

subplot(1,2,2)
plot(z,Sgen)
xlabel ('x [m]')
title('Entropia generada')
ylabel('s_{gen} [w/m^3 K]')
grid on
xlim([0,z(end)])

end

