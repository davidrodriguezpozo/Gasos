%% Gràfiques de convergència
clear all;
 

N = [10;50;100;500;750;1000;2500;5000;7500;10000];

N_2 = [100;500;750;1000;2500;5000;7500;10000];

%Mach in = 0.7

T_1 = [391.21;390.2019;390.0711;390.0237;390.0226;390.0221957;390.0217648;390.0217032;390.0217;390.021687];

ast1 = [390.015;390.015;390.015;390.015;390.015;390.015;390.015;390.015;390.015;390.015];
    
P_1 = [4.5202;4.609;4.619;4.6226;4.6226;4.622678013;4.622710251;4.622714861;4.622715714;4.622716014];

asp1 = [4.6227;4.6227;4.6227;4.6227;4.6227;4.6227;4.6227;4.6227;4.6227;4.6227];

v_1 = [303.594;296.97;296.2361;295.9724;295.9660787;295.9638673;295.9614764;295.9611345;295.9610711;295.961049];

asv1 = [295.96;295.96;295.96;295.96;295.96;295.96;295.96;295.96;295.96;295.96];

rho_1 = [4.0259;4.1157;4.1259;4.1296;4.129706427;4.129737284;4.129770647;4.129775417;4.129776301;4.12977661];

asrho1 = [4.1298;4.1298;4.1298;4.1298;4.1298;4.1298;4.1298;4.1298;4.1298;4.1298];

Mach_1 = [0.76467;0.7498;0.7482;0.7476;0.747625161;0.747623476;0.747624186;0.747625494;0.74762605;0.747626349];

asmach1 = [0.7476;0.7476;0.7476;0.7476;0.7476;0.7476;0.7476;0.7476;0.7476;0.7476];

%Machin = 0.78779

T_2 = [371.2090367;364.1641595;363.9598951;376.1338252;376.1398693;376.1406916;376.1408439;376.1408972];

ast2 = [376.14;376.14;376.14;376.14;376.14;376.14;376.14;376.14];
    
P_2 = [3.82570101;3.586823646;3.580339914;4.020673867;4.020957168;4.020996028;4.021003226;4.021005746];

asp2 = [4.021;4.021;4.021;4.021;4.021;4.021;4.021;4.021];

v_2 = [383.0566;400.8137534;401.3143688;369.3167304;369.296644;369.2938822;369.2933707;369.2931916];

asv2 = [369.293;369.293;369.293;369.293;369.293;369.293;369.293;369.293];

rho_2 = [3.59096017;3.431870776;3.427589724;3.72455644;3.724759023;3.724786878;3.724792037;3.724793844];

asrho2 = [3.724;3.724;3.724;3.724;3.724;3.724;3.724;3.724];

Mach_2 = [0.997581688;1.04810891;1.049613687;0.949897445;0.949898679;0.949910663;0.949915875;0.949918708];

asmach2 = [0.95;0.95;0.95;0.95;0.95;0.95;0.95;0.95];


%Mach in = 2;

T_3 = [451.6368554;437.1253012;433.8729077;432.5339781;432.500609;432.4888974;432.4762274;432.4744152;432.4740795;432.473962];

ast3 = [432.472;432.472;432.472;432.472;432.472;432.472;432.472;432.472;432.472;432.472];
    
P_3 = [6.211789063;5.926703758;5.864308624;5.838777387;5.83814223;5.8379193;5.837678194;5.837643705;5.837637316;5.83763508];

asp3 = [5.837;5.837;5.837;5.837;5.837;5.837;5.837;5.837;5.837;5.837];

v_3 = [728.6981143;739.209764;741.5162706;742.4603789;742.4838678;742.4921112;742.5010291;742.5023046;742.5025408;742.5026235];

asv3 = [742.50;742.50;742.50;742.50;742.50;742.50;742.50;742.50;742.50;742.50];

Mach_3 = [1.72111292;1.765403538;1.776637124;1.781100591;1.78118411;1.781207325;1.781217686;1.781212121;1.781209261;1.781207642];

asmach3 = [1.7815;1.7815;1.7815;1.7815;1.7815;1.7815;1.7815;1.7815;1.7815;1.7815];

%Mach  in = 2 amb ones de xoc

T_4 = [658.9173068;653.4329057;651.973552;651.3686647;651.3558145;651.3518093;651.3486935;651.3488281;651.3490035;651.3491157];

ast4 = [651.34;651.34;651.34;651.34;651.34;651.34;651.34;651.34;651.34;651.34];
    
P_4 = [21.0189804;21.43285787;21.52136422;21.56245124;21.56396212;21.56459494;21.56552659;21.56577766;21.5658547;21.56589204];

asp4 = [21.5659;21.5659;21.5659;21.5659;21.5659;21.5659;21.5659;21.5659;21.5659;21.5659];

v_4 = [314.1913613;305.559568;303.6233387;302.7636293;302.7364435;302.7256981;302.7111721;302.7077103;302.70671;302.7062387];

asv4 = [302.706;302.706;302.706;302.706;302.706;302.706;302.706;302.706;302.706;302.706];

Mach_4 = [0.610428298;0.596234818;0.593154733;0.591799777;0.591757222;0.591740423;0.591717752;0.591712363;0.591710807;0.591710074];

asmach4 = [0.5915;0.5915;0.5915;0.5915;0.5915;0.5915;0.5915;0.5915;0.5915;0.5915];













figure;

subplot(2,2,1);
semilogx(N,T_1); hold on; grid on;
semilogx(N,ast1,'--');
title('Convergencia de les temperatures per $M_{in} = 0.7$');
xlabel('Volums de control');
ylabel('Temperatura final [K]');
ylim([389.8 391.5]);


subplot(2,2,2);
semilogx(N,P_1); hold on; grid on;
semilogx(N,asp1,'--');
title('Convergencia de les pressions per $M_{in} = 0.7$');
xlabel('Volums de control');
ylabel('Pressi\''o final [bar]');

subplot(2,2,3);
semilogx(N,v_1); hold on; grid on;
semilogx(N,asv1,'--');
title('Convergencia de les velocitats per $M_{in} = 0.7$');
xlabel('Volums de control');
ylabel('Velocitat final [m/s]');
ylim([295.5 304]);


subplot(2,2,4);
semilogx(N,Mach_1); hold on; grid on;
semilogx(N,asmach1,'--');
title('Convergencia del Mach per $M_{in} = 0.7$');
xlabel('Volums de control');
ylabel('Mach final');


figure;

subplot(2,2,1);
semilogx(N_2,T_2); hold on; grid on;
semilogx(N_2,ast2,'--');
title('Convergencia de les temperatures per $M_{in} = 0.78799$');
xlabel('Volums de control');
ylabel('Temperatura final [K]');
ylim([360 378]);

subplot(2,2,2);
semilogx(N_2,P_2); hold on; grid on;
semilogx(N_2,asp2,'--');
title('Convergencia de les pressions per $M_{in} = 0.78799$');
xlabel('Volums de control');
ylabel('Pressi\''o final [bar]');

subplot(2,2,3);
semilogx(N_2,v_2); hold on; grid on;
semilogx(N_2,asv2,'--');
title('Convergencia de les velocitats per $M_{in} = 0.78799$');
xlabel('Volums de control');
ylabel('Velocitat final [m/s]');
ylim([365 405]);

subplot(2,2,4);
semilogx(N_2,Mach_2); hold on; grid on;
semilogx(N_2,asmach2,'--');
title('Convergencia del Mach per $M_{in} = 0.78799$');
xlabel('Volums de control');
ylabel('Mach final');
ylim([0.9 1.1]);


figure;

subplot(2,2,1);
semilogx(N,T_3); hold on; grid on;
semilogx(N,ast3,'--');
title('Convergencia de les temperatures per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Temperatura final [K]');

subplot(2,2,2);
semilogx(N,P_3); hold on; grid on;
semilogx(N,asp3,'--');
title('Convergencia de les pressions per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Pressi\''o final [bar]');

subplot(2,2,3);
semilogx(N,v_3); hold on; grid on;
semilogx(N,asv3,'--');
title('Convergencia de les velocitats per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Velocitat final [m/s]');

subplot(2,2,4);
semilogx(N,Mach_3); hold on; grid on;
semilogx(N,asmach3,'--');
title('Convergencia del Mach per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Mach final');
ylim([1.72 1.79]);

figure;

subplot(2,2,1);
semilogx(N,T_4); hold on; grid on;
semilogx(N,ast4,'--');
title('Convergencia de les temperatures per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Temperatura final [K]');


subplot(2,2,2);
semilogx(N,P_4); hold on; grid on;
semilogx(N,asp4,'--');
title('Convergencia de les pressions per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Pressi\''o final [bar]');

subplot(2,2,3);
semilogx(N,v_4); hold on; grid on;
semilogx(N,asv4,'--');
title('Convergencia de les velocitats per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Velocitat final [m/s]');

subplot(2,2,4);
semilogx(N,Mach_4); hold on; grid on;
semilogx(N,asmach4,'--');
title('Convergencia del Mach per $M_{in} = 2$');
xlabel('Volums de control');
ylabel('Mach final');
