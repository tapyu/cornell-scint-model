function [PeVec] = estimatePe(S4,tau0Vec,C_N0)
% estimatePe  Estimate the DPSK bit error probability based on the S4 and
%             tau0 of the continuous-time baseband scintillation signal
%             z(t), where z(t) is assumed to be normalized such that
%             E[|z(t)|^2] = 1.
%
% INPUTS
% S4           S4 index of the continuous-time baseband scintillation 
%              signal z(t).
%
% tau0Vec      Nm-by-1 vector of values of tau0 of the continuous-time baseband
%              scintillation signal z(t), in seconds.  Each element of
%              tau0Vec corresponds to the value of tau0 for one of Nm
%              different candidate models.  These are:
%              (1) Gaussian model
%              (2) 2nd-order Butterworth model
%              (3) f4 model
%
% C_N0         Carrier-to-noise ratio of the received signal, in dB-Hz
%
% OUTPUTS
% PeVec        Nm-by-1 vector of bit error probabilities.
% 
%+------------------------------------------------------------------+
% References:  Based on the closed-form solution for fast Ricean fading
%              in Simon and Alouini Sec. 8.2.5.2.  See also my papers on
%              scintillation characterization and modeling (~2007).
%
% Author:  Todd Humphreys, July 31, 2007
%+==================================================================+

%----- Preliminaries
cnr = 10^(C_N0/10);
Tb = 0.02;
Eb_N0 = cnr*Tb;
Nm = length(tau0Vec);
if(Nm ~= 3) error('Nm must equal 3'); end
PeVec = zeros(Nm,1);

%----- Solve for the K value corresponding to the continuous-time z(t)
if(S4>0)
  m = max(1,1/(S4^2)); 
  Kprime = sqrt(m^2 - m)/(m - sqrt(m^2 - m));
else
  PeVec = zeros(Nm,1);
  return;
end

%----- Solve for sigmaxi2 and zbar
sigmaxi2 =  1/(2*(1 + Kprime));
zbar = sqrt(Kprime*2*sigmaxi2);

%----- Evaluate sigma2 and rho
[sigma2Vec,rhoVec] = evalSigmaAndRho(tau0Vec,sigmaxi2);

%----- Evaluate K and Omega
KVec = zbar^2./(2*sigma2Vec);
OmegaVec = zbar^2 + 2*sigma2Vec;
gammaBarVec = OmegaVec*Eb_N0;

%----- Evaluate Pe for all models
PeVec = 0.5*((1 + KVec + gammaBarVec.*(1 - rhoVec))...
             ./(1 + KVec + gammaBarVec))...
        .*exp(-KVec.*gammaBarVec./(1 + KVec + gammaBarVec));

