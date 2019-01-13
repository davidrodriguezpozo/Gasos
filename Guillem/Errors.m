function [ ] = Errors(N,Err_mom, Err_mass,Err_mom)


figure()
plot(N,Err_mom)
grid on
xlabel('Número de volums de control')
ylabel('Error momentum')
legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')

figure()
plot(N,Err_mass)
grid on
xlabel('Número de volums de control')
ylabel('Error massa')
legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')

figure()
plot(N,Err_ene)
grid on
xlabel('Número de volums de control')
ylabel('Error energia')
legend('M=0.1','M=0.2','M=0.3','M=0.4','M=0.5','M=0.6','M=0.7')

end

