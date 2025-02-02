function Deltap = Body_Water_Model(sample,theta,xdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Triple Oxygen Isotope Animal Body Water             %
%       By Benjamin H. Passey, modified and coded by Huanting Hu      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%% Herbivore (Emu) %%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Model Input Parameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Animal-Specific Parameters
Mass=84;                    % Animal Mass (kg) 40 kg for emus
Met_pre_exp=2.955;          % Metabolic pre-exponent
Met_exp=0.727;              % Metabolic exponent
O2_conv_fact=0.00216;       % Oxygen conversion factor (mol/KJ)
WEI=theta(5);                   % Water economy index (ml/KJ) for ostrich (Williams et al., 1993)
FecalH2O=0.6;               % Fecal H2O content (%)
Sweat_Vapor_ratio=0.5;      % Sweat/(Sweat+Vapor)
Frac_used_O2=0.2;           % Fraction of O2 used
Z_value=10.5;               % Z-value
Body_T=38+273.15;           % Body temperature (K)

% Environmental Parameters
rh = theta(3);
T=273.15+15;                % Enviroment temperature (K)
d18O_mw=theta(2);                 % Meteoric water d18O composition (per mil) Namib
D17O_mw=0.09;               % Meteoric water D17O composition (per mil)
d18O_atmO2=theta(6);          % Atmospheric O2 d18O composition (23.881, Barkan and Luz, 2011)
D17O_atmO2=theta(1);          % Atmospheric O2 D17O composition (-0.506, Barkan and Luz, 2011)

% Food Parameters
Digest=0.7;                 % Relative Digestibility
Energy_effi=0.9;            % Energy extraction efficiency
Carb_content=0.85;          % Food carbon content
Fat_content=0.05;           % Food fat content
Protein_content=0.10;       % Food protein content
H2Oassoc_content=0.55;      % Food associated H2O content (Williams et al., 1993)
Plant_Meat_ratio=1.00;      % Plant/(Plant+Meat) for different diet
Leaf_ratio=theta(4);            % Ingested leaf ratio relative to total ingested vegetation

%%%%%%%%%%%%%%%% Do not Change Anything after this line %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Constants         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lamda=0.528;                % Reference lamda value for water system
R_18_SMOW=0.0020052;        % O-18 isotope ratio, standard SMOW
R_17_SMOW=0.00038;          % O-17 isotope ratio, standard SMOW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Leaf Water Model                   %
%          By J. Roden, modified by C. Tai          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs for leaf water
a_vap_mw=1/alpha1(T);       % Fractionation factor between water vapor and meteoric water (Function:alpha1)
d18O_atm_vap=a_vap_mw*(d18O_mw+1000)-1000;   
                            % Atomospheric humidity, d18O
RH=rh*100;                  % Relative humidity (%)
T_diff=0;                   % Difference between leaf and air temperature (K)
baro_press=95;              % Barometric pressure (kPa) Elevation fuction, can be obtained from talbes
stoma_cond=0.3;             % Stomatal conductance range from 0.001 to 1 mol/m2/s
boud_cond=1;                % Boundary layer conductance range from 0.2 to 3 mol/m2/s

% Constants for models
alpha_k=1.032;              % Kinetic fracionation (Cappa et al., 2003)
alpha_kb=1.021;             % Kinetic fractionation at boundary layer (Cappa et al., 2003)
eO_auto=27;                 % Autotrophic fractionation (sucrose) (Yakir and DeNiro,1990; Luo and Sternberg, 1992; Sternberg and DeNiro, 1983)
eO_hetero=27;               % Heterotrophic fractionation (Cellulos) (Yakir and DeNiro,1990; Luo and Sternberg, 1992; Sternberg and DeNiro, 1983)
fO_ex_mew=0.42;             % Fraction of exchange with medium water (Roden and Ehleringer, 1999)

% Calculations
Leaf_T=T+T_diff;            % Leaf temperature calculated
e_sat=(101325*exp((((-0.1299*(1-(373.15/T))-0.6445)*(1-(373.15/T))-1.976)*(1-(373.15/T))+13.3185)*(1-(373.15/T))))/1000;
                            % Saturation vapor pressure (kPa)
e_a=rh*e_sat;               % Ambient vapor pressure (kPa)
e_i=(101325*exp((((-0.1299*(1-(373.15/Leaf_T))-0.6445)*(1-(373.15/Leaf_T))-1.976)*(1-(373.15/Leaf_T))+13.3185)*(1-(373.15/Leaf_T))))/1000;
                            % Leaf vapor pressure (kPa)
g=1/(1/stoma_cond+1/boud_cond); 
                            % Total leaf conductance to water vapor (mol/m2/s)
E=((e_i-e_a)/baro_press)*g; % Leaf transpiration (mol/m2/s)
W_i=e_i/baro_press;         % Water vapor mole fraction
W_a=e_a/baro_press;         % Ambient water vapor mole fraction
W_s=((stoma_cond*W_i)-E*(1-W_i/2))/(stoma_cond-E/2);
                            % Leaf surface water vapor, mole fraction Ball (1987) in Stomatal Function pp445-476
e_s=W_s*baro_press;         % Vapor pressure at leaf surface (kPa)
alpha_equi=alpha1(T);       % Equilibrium fractionation (liquid-vapor) at leaf temperature (Function alpha1) 
R_source=R_18_SMOW*(1+d18O_mw/1000);  
                            % O isotope ratio of source water
R_atm_vap=R_18_SMOW*(1+d18O_atm_vap/1000);
                            % O isotope ratio of atmospheric water vapor
R_Leaf=alpha_equi*(alpha_k*R_source*(e_i-e_s)/e_i+alpha_kb*R_source*(e_s-e_a)/e_i+R_atm_vap*e_a/e_i);
                            % O isotope ratio of leaf water
                          
%%%%%%%%%%%%%%%%% 
%     Results   % 
%%%%%%%%%%%%%%%%%

d18O_Leaf=(R_Leaf/R_18_SMOW-1)*1000;
                            % Modeled leaf water d18O compositions (per mil)
d18O_medw=0.5*(d18O_mw+d18O_Leaf);
                            % Modeled medium water for cellulose synthesis (per mil) 
d18O_Cellulos=fO_ex_mew*(d18O_medw+eO_hetero)+(1-fO_ex_mew)*(d18O_Leaf+eO_auto);
                            % Modeled cellulos d18O compositions (per mil)


                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               18O Body water Model                %
%                   By Kohn, 1996                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%     Input    %
%%%%%%%%%%%%%%%%

Energy=(10^Met_pre_exp)*(Mass^Met_exp);
                            % Animal Total Metabolic Energy kJ
                            

% Ingested Food 
Food_energy=(Carb_content*17.3+Fat_content*39.7+Protein_content*20.1)*Digest*Energy_effi*1000;
                            % Food metabolic energy kJ/kg
Food_O2=Carb_content*15.4+Fat_content*2+Protein_content*3;
                            % Moles of O2 from per kg food (mole/kg)
Food_consumed=Energy/Food_energy;
                            % Total consumed food (kg) (different from actual ingested food)
Food_O2_content=Food_consumed*Food_O2*Digest*Energy_effi;                            
                            % Food O2 content ingested (mole) 
Food_H2=Carb_content*30.9+Fat_content*60+Protein_content*11;
                            % Moles of H2 from per kg food (mole/kg)
Food_H2_content=Food_consumed*Food_H2*Digest*Energy_effi;
                            % Food H2 content ingested (mole)
Food_H2_O2=Food_H2_content/2-Food_O2_content;                            
                            % Food H2/2- Food O2 (mole)
Food_H2O_assoc=Food_consumed*55.56/(1-H2Oassoc_content)*H2Oassoc_content;
                            % Influx of unbound food H2O (associated H2O,like leaf and stem water, mole) [Food_consumed*55.6 is the mole of ingested free food water]
Food_H2O_assoc_Leaf=Food_H2O_assoc*Leaf_ratio; 
                            % Mole of leaf water ingested
Food_H2O_assoc_Stem=Food_H2O_assoc-Food_H2O_assoc_Leaf;
                            % Mole of stem water ingested
Food_urea=6*Food_consumed*Protein_content*Digest*Energy_effi;  
                            % Mole of Urea or Uric acid (mole)
%%% Missing from line110 to 118
                            
% Atmospheric O2
O2_resp=Energy*O2_conv_fact;
                            % O2 respired (mol)

% Atmospheric H2O
Air_ex=O2_resp/Frac_used_O2/0.21*22.4;
                            % Amount of air fluxed through the lungs (L)
H2O_sat_content=10^(0.686+0.027*(T-273.15))/760/22.4;
                            % Saturation concentration of H2O in air at ambient temperature (mole/L)
Air_H2O=rh*H2O_sat_content*Air_ex;
                            % Amount of H2O from atmosphere (mole)

% Drinking H2O
WEI_turnover=Energy*WEI/18;
                            % WEI predicted water turnover (mol)
H2O_in_nodrink=Food_H2_content+Air_H2O+Food_H2O_assoc;  
                            % Water intake except for drinking water
                            % Different with Konh's model, correct Food_O2_content to Food_O2_content/2
                            
          %%% Different with Kohn's Model %%%
                            % Determine the amount of drinking water (revise negative values)
if H2O_in_nodrink > WEI_turnover
    H2O_turnover=H2O_in_nodrink;
    r=1;                 % r=1, WEI not determine water turnover
else
    H2O_turnover=WEI_turnover;
    r=0;                 % r=0, WEI determine water turnover
end

if WEI_turnover-Air_H2O-Food_H2O_assoc-Food_H2_content>0
    Drinking=WEI_turnover-Air_H2O-Food_H2O_assoc-Food_H2_content;
    t=1;                 % t=1, Drink water  **** note here change Food_H2_content into Food_O2_content ****
else
    Drinking=0;
    t=0;                 % t=0, Do not drink water
end

%%%%%%%%%%%%%%%%
%    Output    %
%%%%%%%%%%%%%%%%

Breath_H2O=Air_ex*10^(0.686+0.027*38)/760/22.4;
                            % Exhaled water vapor from breathing (mole)
Oral_H2O_breath=Breath_H2O/2;      
                            % Oral water loss through breathing                           
M_Fecal=Food_consumed*(1-Digest);
                            % Dry fecal output (kg)
Fecal_O2=0;                 % Fecal O2 content (mole)
Fecal_H2O=M_Fecal/(1-FecalH2O)*FecalH2O*55.56;
                            % Fecal H2O content calculated from dry fecal output (mole)
Urea_O2=Food_urea/2;        % Urea/ureic acid O2 content (mole)
Urine=0.25*H2O_turnover;    % Urinary water lost based on total water turnover (mole)                            
Nasal_H2O=Breath_H2O/4;     % Nasal water vapor loss (mole)
Skin_H2O=1.44*Mass^(2/3);   % Skin water loss (mole)
Residual_H2O=H2O_turnover-Fecal_H2O-Urine-Nasal_H2O-Skin_H2O-Oral_H2O_breath;
                            % Amount of water used for heat loss
Oral_H2O=Residual_H2O*(1-Sweat_Vapor_ratio)+Oral_H2O_breath;    
                            % Total oral water loss
Sweat_H2O=Residual_H2O*Sweat_Vapor_ratio;
                            % Water loss from sweat (liquid)
CO2=O2_resp-Food_H2_O2-Urea_O2;
                            % CO2 loss, corrected for the H2_O2 difference and Urea or uric acid (mole)
RQ=CO2/O2_resp;             % CO2/O2 ratio  
%%% Line80-87 not recorded 


%%%%%%%%%%%%%%%%
%   Fractions  %
%%%%%%%%%%%%%%%%

% Total
T_in=Food_O2_content*2+Food_H2O_assoc_Leaf+Food_H2O_assoc_Stem+O2_resp*2+Air_H2O+Drinking;
                            % Total oxygen input (mole)
T_out_pre=Fecal_H2O+Urea_O2*2+Urine+Nasal_H2O+Skin_H2O+Oral_H2O+Sweat_H2O+CO2*2;
                            % Total oxygen output (mole)
Sec_resi=T_in-T_out_pre;    % Output residual (mole)
T_out=T_out_pre+Sec_resi;   % Total oxygen output revised (mole)
                            
% Input fractions
ffo=Food_O2_content*2/T_in;         % Fraction of food O2 content
flw=Food_H2O_assoc_Leaf/T_in;       % Fraction of leaf water content
fstw=Food_H2O_assoc_Stem/T_in;      % Fraction of stem water content
fO2=O2_resp*2/T_in;                 % Fraction of respired O2 content
fvap=Air_H2O/T_in;                  % Fraction of air H2O
fdw=Drinking/T_in;                  % Fraction of drinking H2O
F_in=ffo+flw+fstw+fO2+fvap+fdw;     % Total input flux (=1)

% Output fractions
ffecO=Fecal_O2/T_out;               % Fraction of fecal O2
ffec=Fecal_H2O/T_out;               % Fraction of fecal H2O
furO=2*Urea_O2/T_out;               % Fraction of urea/uric O2
fur=Urine/T_out;                    % Fraction of urine H2O
fnas=Nasal_H2O/T_out;               % Fraction of nasal H2O
ftrans=Skin_H2O/T_out;              % Fraction of skin transcutaneous H2O
foral=Oral_H2O/T_out;               % Fraction of oral loss H2O
fsw=Sweat_H2O/T_out;                % Fraction of sweat H2O
fresnf=Sec_resi/T_out;              % Fraction of second residul O content
fCO2=CO2*2/T_out;                   % Fraction of exhaled CO2
F_out=ffecO+ffec+furO+fur+fnas+ftrans+foral+fsw+fresnf+fCO2;
                                    % Total output flux (=1)

%%% Line 123-141 and D65 not recorded


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            17O-enabled Body water Model           %
%            By Ben Passey & Huanting Hu            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% Fractionation %
%    Factors    %
%%%%%%%%%%%%%%%%%

%%% Input %%% 
% Fractionation factors bewtween atmpospheric H2O vapor and mw
a_vap_mw_18=1/alpha1(T);                            % alpha_18O
l_vap_mw=0.529;                                     % fractionation exponent between water vapor and liquid water, Barkan and Luz, 2005
a_vap_mw_17=exp(l_vap_mw*reallog(a_vap_mw_18));     % alpha_17O

% Fractionation factors between drinking water and mw                            
a_dw_mw_18=1;                                       % alpha_18O
a_dw_mw_17=1;                                       % alpha_17O

% Fractionation factors between leaf water and mw                            
a_lw_mw_18=(1000+d18O_Leaf)/(1000+d18O_mw);         % alpha_18O
l_lw_mw=-0.0078*rh+0.5216;                          % Fractionation exponent between leaf H2O and meoric H2O, Landais et al, 2006
a_lw_mw_17=exp(l_lw_mw*reallog(a_lw_mw_18));        % alpha_17O

% Fractionation factors between stem water and mw                            
a_stw_mw_18=1;                                      % alpha_18O
a_stw_mw_17=1;                                      % alpha_17O

% Fractionation factors between food cellulos and mw 
a_cell_mw_18=(1000+d18O_Cellulos)/(1000+d18O_mw);   % alpha_18O
l_cell_mw=0.5275;                                   % Fractionation exponent bewteen cellulose and meteoric water, Pack et al, 2012
a_cell_mw_17=exp(l_cell_mw*reallog(a_cell_mw_18));  % alpha_17O

% Fractionation factors for atmospheric O2
a_in_O2_18=(1000+(d18O_atmO2-Z_value)/(1-Frac_used_O2))/(1000+d18O_atmO2);   % alpha_18O
l_in_O2=0.5179;                                     % Fractionation exponent between lung-absorbed O2 and atm.O2, Barkan and Luz, 2005
a_in_O2_17=exp(l_in_O2*reallog(a_in_O2_18));        % alpha_17O

%%% Output %%% 
% Fractionation factors between fecal water and bw 
a_fec_bw_18=1;                                      % alpha_18O
a_fec_bw_17=1;                                      % alpha_17O

% Fractionation factors between urinary water and bw 
a_ur_bw_18=1;                                       % alpha_18O
a_ur_bw_17=1;                                       % alpha_17O

% Fractionation factors between sweat water and bw 
a_sw_bw_18=1;                                       % alpha_18O
a_sw_bw_17=1;                                       % alpha_17O

% Fractionation factors between nasal water vapor and bw
a_nas_bw_18=1/alpha1((T+Body_T)/2);                 % alpha_18O
l_nas_bw=0.529;                                     % Fractionation exponent bewteen nasal water vapor and body water, Barkan and Luz, 2005
a_nas_bw_17=exp(l_nas_bw*reallog(a_nas_bw_18));     % alpha_17O

% Fractionation factors between transcutaneously-lost water vapor and bw
a_trans_bw_18=1/1.018;                              % alpha_18O, Kohn (1996), pg. 1484
l_trans_bw=0.5235;                                  % Estimate is intermediate between binary diffusion of water vapor through atmospheric gas (0.518) and equilibrium vapor - liquid fractionation (0.529).
a_trans_bw_17=exp(l_trans_bw*reallog(a_trans_bw_18)); % alpha_17O

% Fractionation factors between exhaled CO2 and bw
a_CO2_bw_18=alpha2(Body_T);                         % alpha_18O, Fuction alpha2, O'Neil & Adami (1969)
l_CO2_bw=0.5248;                                     % Fractionation exponent between CO2 and H2O, Cao and Liu (2010)
a_CO2_bw_17=exp(l_CO2_bw*reallog(a_CO2_bw_18));     % alpha_17O

% Fractionation factors between exhaled H2O and bw
a_brea_bw_18=1/alpha1(Body_T);                      % alpha_18O
l_brea_bw=0.529;                                    % Fractionation exponent between breath water vapor and body water, Barkan and Luz (2005)
a_brea_bw_17=exp(l_brea_bw*reallog(a_brea_bw_18));  % alpha_17O

%%%%%%%%%%%%%%%%
% Calculations %
%%%%%%%%%%%%%%%%

% Air O2 composition
d17O_atmO2=(exp((D17O_atmO2+lamda*reallog(d18O_atmO2/1000+1)*1000)/1000)-1)*1000;
                            % d17O of atmospheric O2 based on input information
R_18_atmO2=((d18O_atmO2/1000)+1)*R_18_SMOW;
                            % Calculated R_18O of atmospheric O2
R_17_atmO2=((d17O_atmO2/1000)+1)*R_17_SMOW;
                            % Calculated R_17O of atmospheric O2

% Meteoric H2O composition
d17O_mw=(exp((D17O_mw+lamda*reallog(d18O_mw/1000+1)*1000)/1000)-1)*1000;
                            % d17O of meteoric water based on input information
R_18_mw=((d18O_mw/1000)+1)*R_18_SMOW;
                            % Calculated R_18O of meteoric H2O
R_17_mw=((d17O_mw/1000)+1)*R_17_SMOW;
                            % Calculated R_17O of meteoric H2O

%%% Mass Balance Equation %%%

%%% 18R %%%
% Input
R_18_in=R_18_mw*(fvap*a_vap_mw_18+fdw*a_dw_mw_18+flw*a_lw_mw_18+fstw*a_stw_mw_18+ffo*a_cell_mw_18);
                            % Oxygen_18 input term exclude atmospheric O2
R_18_inatm=R_18_atmO2*fO2*a_in_O2_18;
                            % Oxygen_18 input term from atmospheric O2
% Output
r_18_out=ffec*a_fec_bw_18+furO*a_ur_bw_18+fur*a_ur_bw_18+fnas*a_nas_bw_18+ftrans*a_trans_bw_18+foral*a_brea_bw_18+fsw*a_sw_bw_18+fCO2*a_CO2_bw_18+fresnf*1;
                            % Oxygen_18 Output term devided by R_18_bw (Or the pre-constant C in equation: A*R_mw+B*R_atmO2=C*R_bw) 

%%% 17R %%%
% Input
R_17_in=R_17_mw*(fvap*a_vap_mw_17+fdw*a_dw_mw_17+flw*a_lw_mw_17+fstw*a_stw_mw_17+ffo*a_cell_mw_17);
                            % Oxygen_17 input term exclude atmospheric O2
R_17_inatm=R_17_atmO2*fO2*a_in_O2_17;
                            % Oxygen_17 input term from atmospheic O2
% Output
r_17_out=ffec*a_fec_bw_17+furO*a_ur_bw_17+fur*a_ur_bw_17+fnas*a_nas_bw_17+ftrans*a_trans_bw_17+foral*a_brea_bw_17+fsw*a_sw_bw_17+fCO2*a_CO2_bw_17+fresnf*1;
                            % Oxygen_17 Output term devided by R_17_bw (Or the pre-constant C in equation: A*R_mw+B*R_atmO2=C*R_bw) 


%%%%%%%%%%%%%%%%
%    Results   %
%%%%%%%%%%%%%%%%

R_18_bw=(R_18_in+R_18_inatm)/r_18_out;
                            % Oxygen_18 mass balance
R_17_bw=(R_17_in+R_17_inatm)/r_17_out;
                            % Oxygen_17 mass balance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d18O_bw=(R_18_bw/R_18_SMOW-1)*1000;
                            % d18O of body water
d17O_bw=(R_17_bw/R_17_SMOW-1)*1000;
                            % d17O of body water
deltap  = xdata(sample,2);
Deltap=(reallog(d17O_bw/1000+1)-lamda*reallog(deltap/1000+1))*1000;
                            % D17O of body water

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equivalent carbonates information
%d18O_carb_SMOW=1.027*(d18O_bw+1000)-1000;       % Carbonate d18O in SMOW
%d18O_CO2_SMOW=1.008*(d18O_carb_SMOW+1000)-1000; % Acid digested CO2 d18O
%d18O_carb_PDB=(d18O_carb_SMOW-30.86)/1.03086;   % carbonate d18O in PDB

% relative humidity range
%end

%%%%%%%%%%%%%%%%
%     Plots    %
%%%%%%%%%%%%%%%%

% X axis can be changed into different parameters
%scatter(d18O_bw,D17O_bw,20,Rel_Humidity,'filled');
%colorbar;
%box on;
% WEI=linspace(0.01,1,100);  % X-axis
% figure
% subplot(2,1,1) 
% plot(WEI,d18O_bw,'b')
% xlabel('WEI');
% ylabel('d18O_b_w');
% axis([0 1 -10 10])
% subplot(2,1,2)
% figure
% plot(WEI,D17O_bw,'b')
% xlabel('WEI');
% ylabel('D17O_b_w');
% axis([0 1 -0.2 0.1])

%%%%%%%%%%%%%%%
%    File     %
%%%%%%%%%%%%%%%

%j=length(d18O_bw);
%File=nan(j,5);
%File(:,1)=d18O_bw;
%File(:,2)=D17O_bw;
%File(:,3)=WEI;
%File(:,4)=r;
%File(:,5)=t;
%save('WEI_new.mat','File');
%dlmwrite('WEI_new.csv',File);
