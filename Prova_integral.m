
Tin = 400;
Trs = 300;

fun = @(var) 1022 - 0.166*var + 3.5025e-4*var.^2;

Cpr = integral(fun,Tin,Trs)/(Trs-Tin);