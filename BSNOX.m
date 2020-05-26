function [ w_NO ] = BSNOX( T_NO,P_atm,lambda,P_NO,T_BDC,P_BDC,P_EXH)
R_u=8315;           % Univerzální plznová konstanta
psi=3.773;          % Molární poměr N/O 
y=18/8;             % Molární poměr H/C (Použit Iso-Octane)
epsilon=4/(4+y);    % y je molární poměr H/C

% Vypočítat rovnovážnou konstantu při dané teplotě (posun vodního plynu)
K_wgs=exp(2.743-1.761*10^3/T_NO-1.611*10^6/(T_NO^2)+.2803*10^9/(T_NO^3));

% Atomová rovnováha založená na nadměrném množství vzduchu

if 1/lambda<1
    n_CO2=epsilon*(1/lambda);
    n_H2O=2*(1-epsilon)*(1/lambda);
    n_CO=0;
    n_H2=0;
    n_O2=1-(1/lambda);
    n_N2=psi;
    n_b= (1-epsilon)*(1/lambda)+1+psi;
end

if 1/lambda>=1
    A=(K_wgs-1);
    B=-K_wgs*(2*((1/lambda)-1)+epsilon*(1/lambda))+2*(1-epsilon*(1/lambda));
    C=2*K_wgs*epsilon*(1/lambda)*((1/lambda)-1);
    % Sledujte kvadratickou rovnici, moly musí být pozitivní!
    c=(-B-sqrt(B^2-4*A*C))/(2*A);
    n_CO2=epsilon*(1/lambda)-c;
    n_H2O=2*(1-epsilon*(1/lambda))+c;
    n_CO=c;
    n_H2=2*((1/lambda)-1)-c;
    n_O2=0;
    n_N2=psi;
    n_b = (2-epsilon)*(1/lambda)+psi;
end
%__________________________________________________________________________
% Vypočítejte molární zlomky každého prvku
x_CO2=n_CO2/n_b;    x_H2O=n_H2O/n_b;    x_CO=n_CO/n_b;  x_H2=n_H2/n_b;
x_O2=n_O2/n_b;  x_N2_e=n_N2/n_b; 

if 1/lambda>1
    n_prod=1.5;
    z=(T_NO-2.3*10^3)/(7.6*10^2);
    K_p_CO2 = -.55*z^3+1.5*z^2-3*z+9.1;
    Z=(T_NO-2.3*10^3)/(7.6*10^2);
    K_p_CO = -.15*Z^3+.41*Z^2-.92*Z+7.1;
    K_P=10^(K_p_CO2-K_p_CO);
    P_p = (((P_EXH/101325)/n_prod)*(T_NO/T_BDC))^(-1); 
    ALPHA= (2*P_p)/(3*K_P^2) + (((P_p/K_P^2 - (4*P_p^2)/(3*K_P^4) + ...
        (8*P_p^3)/(27*K_P^6))^2 +((4*P_p)/(3*K_P^2) -...
        (4*P_p^2)/(9*K_P^4))^3)^(1/2) + P_p/K_P^2 - (4*P_p^2)/(3*K_P^4) +...
        (8*P_p^3)/(27*K_P^6))^(1/3) - ((4*P_p)/(3*K_P^2) - ...
        (4*P_p^2)/(9*K_P^4))/(((P_p/K_P^2 -(4*P_p^2)/(3*K_P^4) + ...
        (8*P_p^3)/(27*K_P^6))^2 + ((4*P_p)/(3*K_P^2) - ...
        (4*P_p^2)/(9*K_P^4))^3)^(1/2) + P_p/K_P^2 - (4*P_p^2)/(3*K_P^4) + ...
        (8*P_p^3)/(27*K_P^6))^(1/3);
    x_O2=(ALPHA/(2*n_prod))*(x_CO2*((1-ALPHA)/n_prod)); 
    %x_CO=x_CO+(ALPHA/n_prod);
end
% Vypočítat rovnovážné koncentrace
X_O2_e=x_O2*P_BDC/(R_u*T_NO);
%__________________________________________________________________________
% Rovnovážná konstanta pro reakci kyslíku na kyslík
Kp_7=3.6*10^3*exp(-31090/T_NO)*318.3; %318.3 converts atm^(1/2) to pa^(1/2)
x_O_e= (Kp_7*X_O2_e^(1/2))/((R_u*T_NO)^(1/2))/(P_BDC/(R_u*T_NO)); %kmol/m^3
% Konstantní rychlost reakce vpřed  [m^3/kmol-s]
k_1f=1.82*10^11*exp(-38370/T_NO);
% Vypočítat změnu NO koncentrace jako funkci času
dNOdt=2*k_1f*x_O_e*x_N2_e*P_BDC/(R_u*T_NO);
% Vypočítat dobu pobytu
t_NO=(8*10^(-16)*T_NO*exp(58300/T_NO))/(P_NO/101325)^(1/2);
% Vypočítat NO PPM
PPM_NO = dNOdt*t_NO*10^6;
%__________________________________________________________________________
%Converting to mass fractions of constituents:
m_CO2=x_CO2*(44.01);    m_H2O=x_H2O*(18.015);    m_CO=x_CO*(28.01);
m_H2=x_H2*(2.01);  m_O2=x_O2*32;    m_N2_e=x_N2_e*28;

%Adding mass fractions of all constituents:

m_all=m_CO2+m_H2O+m_CO+m_H2+m_O2+m_N2_e;

%Calculating mass fraction of NO:

w_NO=dNOdt*t_NO*30.01/m_all;

end