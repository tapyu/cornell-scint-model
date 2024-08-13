function [sigma2Vec,rhoVec] = evalSigmaAndRho(tau0Vec,sigmaxi2)
% evalSigmaAndRho   Evaluate the quantities sigma2 and rho used in calculating
%                   the bit error probability of DPSK detection for three
%                   different candidate models of 
%                   Rxi(tau) = (1/2)*E[conj(xi(t))*xi(t+tau)]
%
% INPUTS
% tau0Vec      Nm-by-1 vector of optimal values of tau0.  Each element of
%              tau0Vec corresponds to the optimal for one of Nm different
%              candidate models.  These are:
%              (1) Gaussian model
%              (2) 2nd-order Butterworth model
%              (3) f4 model
%
% sigmaxi2     The value of Rxi(tau) at tau = 0.
%
% OUTPUTS
% sigma2Vec    Nm-by-1 vector of sigma2 values corresponding to the
%              models whose tau0 values are stored in tau0Vec.
%
% rhoVec       Nm-by-1 vector of rho values corresponding to the
%              models whose tau0 values are stored in tau0Vec.
%
% 
%+------------------------------------------------------------------+
% References:  Roger Dana, "Effects of Ionospheric Scintillation on
%              Differential Demodulation of GPS data," IEEE Trans on
%              Aerospace and Electronic Systems, vol. 33 no. 3,
%              Jul. 1997.  
%
%              See also my notes "Computation of sigma and rho for
%              2nd-order Butterworth," from July 31, 2007. 
%
% Notes:       The "brute force" integration techniques (see bottom of
%              this file) agree very well with the theoretical results
%              for sigma2 and rho.
%
% Author:  Todd Humphreys, July 31, 2007
%+==================================================================+

%----- Preliminaries
Nm = length(tau0Vec);
if(Nm ~= 3) error('Nm must equal 3'); end
% dtau defines the resolution used in evaluating the integrals
Tb = 0.02;
rhoVec = zeros(Nm,1);
sigma2Vec = zeros(Nm,1);

%----- Gaussian model
mm = 1;
tau0 = tau0Vec(mm);
betag = tau0/Tb;
I11g = (sqrt(pi)*betag/2)*erf(1/betag);
I12g = (sqrt(pi)*betag/2)*erf(2/betag); 
I21g = (betag^2/2)*(1 - exp(-1^2/betag^2));
I22g = (betag^2/2)*(1 - exp(-2^2/betag^2));
R1g = (2*I11g - 2*I21g);
R1mR2g = (4*I11g - 2*I12g + I22g - 4*I21g);
R2g = -(R1mR2g - R1g);
sigma2Vec(mm) = sigmaxi2*R1g;
rhoVec(mm) = R2g/R1g;

%----- 2nd-order Butterworth model
mm = 2;
beta = 1.23964643681047; % Eq. (7)
tau0 = tau0Vec(mm);
q = (beta*Tb/tau0); % Eq. (7)
fq = exp(-q)*(cos(q)-sin(q));
f2q = exp(-2*q)*(cos(2*q)-sin(2*q));
sigma2Vec(mm) = (sigmaxi2/q^2)*(2*q + fq - 1);
sigma2_rho = (sigmaxi2/(2*q^2))*(f2q - 2*fq + 1);
rhoVec(mm) = sigma2_rho/sigma2Vec(mm);    

%----- f4 model
mm = 3;
a4 = 2.14619322062058;
tau0 = tau0Vec(mm);
betaf4 = tau0/a4/Tb;
I11f4 = 2*betaf4 - (1 + 2*betaf4)*exp(-1/betaf4);
I12f4 = 2*betaf4 - (2 + 2*betaf4)*exp(-2/betaf4);
I21f4 = 3*betaf4^2 - (1^2 + 3*1*betaf4 + 3*betaf4^2)*exp(-1/betaf4);
I22f4 = 3*betaf4^2 - (2^2 + 3*2*betaf4 + 3*betaf4^2)*exp(-2/betaf4);
R1f4 = (2*I11f4 - 2*I21f4);
R1mR2f4 = (4*I11f4 - 2*I12f4 + I22f4 - 4*I21f4);
R2f4 = -(R1mR2f4 - R1f4);
sigma2Vec(mm) = sigmaxi2*R1f4;
rhoVec(mm) = R2f4/R1f4;



% %------------ "Brute Force" Techniques
% dtau = 1e-5;
% tVec = [0:dtau:Tb-dtau]';
% Nt = length(tVec);
% rhoVecBF = zeros(Nm,1);
% sigma2VecBF = zeros(Nm,1);
% %----- Gaussian model
% mm = 1;
% tau0 = tau0Vec(mm);
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - Tb:dtau:t]';
%   Rxi = exp(-(tauVec/tau0).^2);
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% sigma2VecBF(mm) = (1/Tb^2)*S;
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - 2*Tb:dtau:t-Tb]';
%   Rxi = exp(-(tauVec/tau0).^2);
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% rhoVecBF(mm) = (1/Tb^2)*S/sigma2VecBF(mm);
% % Compare Gaussian results with 1-D theoretical integration
% tv = [0:dtau:Tb]';
% y1 = (tau0*sqrt(pi)/(2*Tb^2))*(erf(tv/tau0) + erf((Tb - tv)/tau0));
% y2 = (tau0*sqrt(pi)/(2*Tb^2))*(erf((2*Tb - tv)/tau0) - erf((Tb - tv)/tau0));
% sigma2 = trap_int(y1,dtau);
% rho = trap_int(y2,dtau)/sigma2;
% % Compare Gaussian results with theory
% betag = tau0/Tb;
% I11g = (sqrt(pi)*betag/2)*erf(1/betag);
% I12g = (sqrt(pi)*betag/2)*erf(2/betag); 
% I21g = (betag^2/2)*(1 - exp(-1^2/betag^2));
% I22g = (betag^2/2)*(1 - exp(-2^2/betag^2));
% % R variables
% R1g = (2*I11g - 2*I21g);
% R1mR2g = (4*I11g - 2*I12g + I22g - 4*I21g);
% R2g = -(R1mR2g - R1g);
% sigma2g = R1g;
% rhog = R2g/sigma2;
% % Compare with Dana's Eq. 22
% dxi = 0.0001;
% xiVec = [0:dxi:1]';
% y1Vec = 2*(1 - xiVec).*exp(-(Tb*xiVec/tau0).^2);
% R1eq22 = trap_int(y1Vec,dxi);
% y2Vec = xiVec.*exp(-(Tb*xiVec/tau0).^2);
% xiVec = [1:dxi:2]';
% y3Vec = (2-xiVec).*exp(-(Tb*xiVec/tau0).^2);
% R2eq22 = trap_int(y2Vec,dxi) + trap_int(y3Vec,dxi);

% %----- 2nd-order Butterworth model
% mm = 2;
% beta = 1.23964643681047;
% tau0 = tau0Vec(mm);
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - Tb:dtau:t]';
%   Rxi = exp(-beta*abs(tauVec)/tau0)...
%         .*(cos(beta*tauVec/tau0) + sin(beta*abs(tauVec)/tau0));
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% sigma2VecBF(mm) = (1/Tb^2)*S;
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - 2*Tb:dtau:t-Tb]';
%   Rxi = exp(-beta*abs(tauVec)/tau0)...
%         .*(cos(beta*tauVec/tau0) + sin(beta*abs(tauVec)/tau0));
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% rhoVecBF(mm) = (1/Tb^2)*S/sigma2VecBF(mm);
% % Compare with Dana's Eq. 22
% dxi = 0.00001;
% xiVec = [0:dxi:1]';
% y1Vec = 2*(1 - xiVec).*exp(-beta*abs(Tb*xiVec)/tau0)...
%         .*(cos(beta*Tb*xiVec/tau0) + sin(beta*abs(Tb*xiVec)/tau0));
% R1eq22 = trap_int(y1Vec,dxi);
% y2Vec = xiVec.*exp(-beta*abs(Tb*xiVec)/tau0)...
%         .*(cos(beta*Tb*xiVec/tau0) + sin(beta*abs(Tb*xiVec)/tau0));
% xiVec = [1:dxi:2]';
% y3Vec = (2-xiVec).*exp(-beta*abs(Tb*xiVec)/tau0)...
%         .*(cos(beta*Tb*xiVec/tau0) + sin(beta*abs(Tb*xiVec)/tau0));
% R2eq22 = trap_int(y2Vec,dxi) + trap_int(y3Vec,dxi);
% % Compare with my theoretical integration
% q = (beta*Tb/tau0);
% fq = exp(-q)*(cos(q)-sin(q));
% f2q = exp(-2*q)*(cos(2*q)-sin(2*q));
% sigma2 = (1/q^2)*(2*q + fq - 1);
% sigma2_rho = (1/(2*q^2))*(f2q - 2*fq + 1);
% rho = sigma2_rho/sigma2;    

% %----- f4 model
% mm = 3;
% a4 = 2.14619322062058;
% tau0 = tau0Vec(mm);
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - Tb:dtau:t]';
%   Rxi = (1 + a4*abs(tauVec)/tau0).*exp(-a4*abs(tauVec)/tau0);
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% sigma2VecBF(mm) = (1/Tb^2)*S;
% S = 0;
% for ii=1:Nt
%   t = tVec(ii);
%   tauVec = [t - 2*Tb:dtau:t-Tb]';
%   Rxi = (1 + a4*abs(tauVec)/tau0).*exp(-a4*abs(tauVec)/tau0);
%   S = S + dtau*trap_int(Rxi,dtau);
% end
% rhoVecBF(mm) = (1/Tb^2)*S/sigma2VecBF(mm);
% % Compare f4 results with theory
% betaf4 = tau0/a4/Tb;
% I11f4 = 2*betaf4 - (1 + 2*betaf4)*exp(-1/betaf4);
% I12f4 = 2*betaf4 - (2 + 2*betaf4)*exp(-2/betaf4);
% I21f4 = 3*betaf4^2 - (1^2 + 3*1*betaf4 + 3*betaf4^2)*exp(-1/betaf4);
% I22f4 = 3*betaf4^2 - (2^2 + 3*2*betaf4 + 3*betaf4^2)*exp(-2/betaf4);
% % R variables
% R1f4 = (2*I11f4 - 2*I21f4);
% R1mR2f4 = (4*I11f4 - 2*I12f4 + I22f4 - 4*I21f4);
% R2f4 = -(R1mR2f4 - R1f4);
% sigma2 = R1f4;
% rho = R2f4/sigma2;







