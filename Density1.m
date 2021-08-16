ouput = zeros(10,6);
FracElasticStrain_Final = zeros(10,1);
FracPlasticStrain_Final = zeros(10,1);
for ii =1%:10
% Basic inputs and Model Set-up
% ************************* Model information
Length = 10;
Width = 10;
Shear_G = 25e9;
PoissonRatio = 0.25;
YM = 2*Shear_G*(1+PoissonRatio);
% ************************* Far-field stress tensor
Sigma1_far = 150e6;% Unit: Pa
Sigma3_far = 75e6;
% Pp = (Sigma1_far+Sigma3_far)/2-(Sigma1_far-Sigma3_far)/2/sin(atan(0.6)); % Pore pressure
Pp = 39.6e6;
Theta = 80; % Angle between sigma1 and the damage zone normal vector
fai = (90-Theta)*pi/180;
% ************************* Crack information
k_n = 12000e9; % Pa/m 
k_s = 20e9; % Pa/m
L_s = 4*Shear_G/(pi*(1-PoissonRatio)*k_s); % stiffness length
load(['Realization' num2str(ii) '.mat'])
Crk_num = length(fLeng);
Crk_len = fLeng; % length
Dila_coe = 0.05; % dilatancy factor
Crk_Orien = FracData3(:,5);
rng(1,'twister');
s = rng;
rng(s);
Crk_mu = 0.2*ones(Crk_num,1);
% ************************* Histogram of frictional coefficient
% [~] = Histogram_Statistics( Crk_mu,min(Crk_mu),max(Crk_mu),50 );
%% Stress conditions
Sigma_xx = (Sigma1_far+Sigma3_far)/2+(Sigma1_far-Sigma3_far)*cos(2*fai)/2;
Sigma_yy = (Sigma1_far+Sigma3_far)/2-(Sigma1_far-Sigma3_far)*cos(2*fai)/2;
Sigma_xy = (Sigma1_far-Sigma3_far)*sin(2*fai)/2;
%% Initiation before deviatoric loading
% Initiation of crack information matrix
[Crk_inf,Critical_inf,NonCri_inf] = Crk_Estimate(Crk_num,Crk_len,Crk_mu,Crk_Orien,Shear_G,PoissonRatio,Sigma_xx,Sigma_yy,Sigma_xy,k_s,k_n,Pp,Dila_coe);
%% Matrix strain response
[ Matrix_strain ] = Strain_Solver((Sigma_xx-39.6e6), (Sigma_yy-39.6e6), Sigma_xy, PoissonRatio, YM);
Ini_Strain_xx = Matrix_strain(1,1);
Ini_Strain_yy = Matrix_strain(2,2); % strain_yy will be kept constant
Ini_Strain_xy = Matrix_strain(1,2);
%% stress relaxation -----non-interaction approximation
% Initializing
Xstrain_store = zeros(2,1);
Ystress_store = zeros(2,1);

PR_store = zeros(2,1); % depends on the iterative number
PR_store(1) = PoissonRatio;
YM_store = zeros(2,1);
YM_store(1) = 2*(1+PoissonRatio)*Shear_G;
SM_store = zeros(1,2); % the first colomn is derived from shear stress/strain
SM_store(1,1) = Shear_G;   % the second is derived from YM and PR
SM_store(1,2) = Shear_G;                      
%
% 
% non-critical fractures
[NonCrk_id,~] = size(NonCri_inf);
Epsilon_ElasticCrk = zeros(2,2);
for j = 1:NonCrk_id
    % Extract basic information of each fracture
    Crk_Length = NonCri_inf(j,2);
    alpha = NonCri_inf(j,4);
    Normal_Vector = [cos(alpha);sin(alpha);];
    Disp_UnitVector = NonCri_inf(j,7:8)';
    Mean_disp_Shear = NonCri_inf(j,12);
    Mean_disp_Normal = NonCri_inf(j,13);
    % Boundary strain contributed by each single fracture
    [D_Epsilon_Crk ] = Disp_Strain(Length,Width,Crk_Length,Mean_disp_Shear,Mean_disp_Normal,Normal_Vector,Disp_UnitVector);
    %
    Epsilon_ElasticCrk = Epsilon_ElasticCrk + D_Epsilon_Crk;
    %
end
%
%
% critical fractures
[CriCrk_id,~] = size(Critical_inf);
Delta_Epsilon_PlasticCrk = zeros(2,2);
for i = 1:CriCrk_id
    % Extract basic information of each fracture
    Crk_Length = Critical_inf(i,2);
    alpha = Critical_inf(i,4);
    Normal_Vector = [cos(alpha);sin(alpha);];
    Disp_UnitVector = Critical_inf(i,7:8)';
    Mean_disp_Shear = Critical_inf(i,12);
    Mean_disp_Normal = Critical_inf(i,13);
    % Boundary strain contributed by each single fracture
    [D_Epsilon_Crk ] = Disp_Strain(Length,Width,Crk_Length,Mean_disp_Shear,Mean_disp_Normal,Normal_Vector,Disp_UnitVector);
    %
    Delta_Epsilon_PlasticCrk = Delta_Epsilon_PlasticCrk + D_Epsilon_Crk;
    %
end
%
% calculate strain components
Strain_Crk_input = Epsilon_ElasticCrk+Delta_Epsilon_PlasticCrk; % total fracture strain
YY_strain_temp = Strain_Crk_input(2,2) + Ini_Strain_yy; % total y-strain before relaxation
%
% Solve the increase of sigma_yy by Bisection Method
Delta_Sigmayy_ini = YM/(1-PoissonRatio^2)*Strain_Crk_input(2,2);
[Del_Stress_yy,New_Stress_yy,Mat_Strain_yy,Delta_Epsilon,k,a,b,fp,Crk_inf_final] = Bisection(0,Delta_Sigmayy_ini,Sigma_xx,Sigma_yy,Ini_Strain_yy,Sigma_xy,PoissonRatio,YM,Crk_num,k_s,k_n,Pp,Length,Width,Crk_inf,Strain_Crk_input,Theta,1e-8);
% update stress and strain
New_Strain_xx = Ini_Strain_xx+PoissonRatio*(1+PoissonRatio)/YM*Del_Stress_yy+Delta_Epsilon(1,1);
New_Strain_yy = Ini_Strain_yy;
New_Strain_xy = Ini_Strain_xy+Delta_Epsilon(1,2);
% Update elastic properties
New_Strain=[New_Strain_xx,New_Strain_xy;New_Strain_xy,New_Strain_yy;];
[New_PR, New_YM] = Elasticity_Update(Sigma_xx,New_Stress_yy,New_Strain);
New_SM = Sigma_xy/(2*New_Strain_xy);
New_SM_Cal = New_YM/2/(1+New_PR); % calculated based on Young's modulus and Poisson's ratio
% Store relevant information after one slip (all critical fractures)
Xstrain_store(2) = New_Strain_xx;
Ystress_store(2) = New_Stress_yy;
PR_store(2) = New_PR;
YM_store(2) = New_YM;
SM_store(2,1) = New_SM;
SM_store(2,2) = New_SM_Cal;
%% Final contribution of fracture elastic/plastic deformation
StrainCtrb_ElasticFinal = zeros(2,2);
StrainCtrb_PlasticFinal = zeros(2,2);
for i = 1:Crk_num
    % Extract basic information of each fracture
    Crk_Length = Crk_inf_final(i,2);
    alpha = Crk_inf_final(i,4);
    Normal_Vector = [cos(alpha);sin(alpha);];
    Disp_UnitVector = Crk_inf_final(i,7:8)';
    Mean_disp_Shear = Crk_inf_final(i,12);
    Mean_disp_Normal = Crk_inf_final(i,13);
    % Boundary strain contributed by each single fracture
    [D_Epsilon_Crk ] = Disp_Strain(Length,Width,Crk_Length,Mean_disp_Shear,Mean_disp_Normal,Normal_Vector,Disp_UnitVector);
    %
    if Crk_inf_final(i,6) == 0
        StrainCtrb_ElasticFinal = StrainCtrb_ElasticFinal + D_Epsilon_Crk;
    else
        StrainCtrb_PlasticFinal = StrainCtrb_PlasticFinal + D_Epsilon_Crk;
    end
    %
end
%
FracElasticStrain_Final(ii,1) = StrainCtrb_ElasticFinal(2,2);
FracPlasticStrain_Final(ii,1) = StrainCtrb_PlasticFinal(2,2);
%%
Sigmaxx_final = Sigma_xx;
Sigmayy_final = Ystress_store(end);
Sigmaxy_final = Sigma_xy;

Sigma1_final = (Sigmaxx_final+Sigmayy_final)/2+sqrt((Sigmaxx_final-Sigmayy_final)^2+4*Sigmaxy_final^2)/2;
Sigma3_final = (Sigmaxx_final+Sigmayy_final)/2-sqrt((Sigmaxx_final-Sigmayy_final)^2+4*Sigmaxy_final^2)/2;
PrincipalStress_final = [Sigma1_final Sigma3_final];
theta_final=atan((2*Sigmaxy_final)/(Sigmaxx_final-Sigmayy_final))*180/pi/2;
if theta_final<=0
    theta_final = theta_final+90;
end


% figure
% plot([0 200],[0 0],'k-')
% hold on 
% [~]=Mohr_Circle_line(Sigma1_far/1e6, Sigma3_far/1e6);
% hold on
% [~] = Mohr_Circle_dot(Sigma1_final/1e6, Sigma3_final/1e6);
% hold on
% plot([Pp/1e6 130],[0 (130-Pp/1e6)*0.6],'k-')
% hold on
% axis equal
% ylim([0 50])
% xlim([20 180])
%
ouput(ii,:) = [Sigmayy_final/(1e6),Sigma1_final/(1e6),Sigma3_final/(1e6),(90-theta_final),YM_store(end)/(1e9),PR_store(end),];

end








