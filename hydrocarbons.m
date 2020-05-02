function [ HC ] = hydrocarbons( R_frac,AF_ratio_ac,B,P_peak,imep,C_r,V_d,N_cyl,T_w,N )
%Offset Of Spark Plug From Central Axis of Cylinder
d_splug = 0;
%Calculate Crevice Volume
h_crevice = (3/1000);   %Crevice Height (m)
gap = (1.5/1000);       %Crevice Width
V_crevice = (pi/4)*B^2*h_crevice - (pi/4)*(B-2*gap)^2*h_crevice;
%Calculate Unburned Fraction
f_unburned = (1-R_frac);
%Calculate Fuel Vapor
f_vapor = 1/(1+(AF_ratio_ac/15.09)*14.7);
%Modification Factor Based On Spark Plug
f_mod = (1-.858*(d_splug/B));
%Crevice Emissions Index
SF_crevice = 5443*(P_peak/imep)*(V_crevice/(V_d/N_cyl))*(1/T_w)*...
    f_unburned*f_vapor*f_mod;
%Oil Layer Predictions
P_i = .09875+.00986*imep;
P_ideal = (P_i+P_i*C_r^1.4)/2;
SF_wall = 63024*(1/imep)*(1/((AF_ratio_ac/15.09)*...
    14.7*10^(.0082*T_w)*B))*P_ideal;
%The Threshold of HC Oxidation
P_70 = .209+.0102*imep;
T_70 = 1600+.759*imep-.00051*imep^2;
T_HC = (T_70-T_w)/log(T_70/T_w);
T_HC_adj = 1600+.759*imep-.000051*imep^2;
%Fraction of Cylinder Oxidation
f_ox = 1-(P_70/P_ideal)*(T_HC/T_HC_adj)^3;
RELSP = .829*R_frac/100;
f_ox_ex = .866-.0000146*N-.00007*imep-.007918*RELSP-.0000255*T_w;
HC = (SF_crevice*(1-f_ox)+SF_wall*(1-f_ox))*(f_ox)*(1-f_ox_ex);
end
