function [zScint,zScintA,tScint] = scintModel04(Ts,Nt,tau0,K,Nspa)
% Simulates a realistic complex scintillation time history.
%
% Version 2 --- Outputs both samples spaced by Ts seconds (zScint) and
%               averages over Ts seconds that are also spaced by Ts seconds
%               (zScintA).  Only the underlying high-sample-rate zScintC
%               (effectively continuous time) is normalized so that
%               E[|zScintC(k)|^2] = 1.  Hence, zScint and zScintA are only
%               approximately normalized.
%
% Version 3 --- Like Version 2 but assumes a 2nd-order Butterworth power
%               spectrum given by Sxx(f) = 1/[1 + 16 (f^4/Bd^4)], where Bd =
%               sqrt(2)*beta0/(pi*tau0) is the fading bandwidth, with beta0
%               = 1.23964643681047.  Hence, the scintillaiton model only
%               depends on tau0 and K, nothing else.  Also, the number of
%               samples per average is made an input.
%
% Version 4 --- Like Version 3 but uses an actual 2nd-order Butterworth
%               filter instead of a spectral weighting function to
%               produce the time histories of scintillation.
%
% INPUTS
% Ts           Sampling interval of the output time histories zScint and
%              zScintA, in seconds.
%
% Nt           Number of samples in the output time histories zScint 
%              and zScintA.
%
% tau0         The decorrelation time (in seconds) of the 2nd-order
%              Butterworth model used to generate the scintillation.
%
% K            The Ricean K (noncentrality) parameter
%
% Nspa         Number of sub-samples that are used in calculating the
%              averages in zScintA.  Using more sub-samples increases the
%              accuracy of the averages at the expense of a greater
%              computational burden.
%          
% OUTPUTS
% zScint       Nt-by-1 vector containing the normalized complex
%              scintillation time history in the form of samples with
%              sampling interval Ts.  zScint(k) is the sample taken at
%              time tScint(k).  
%
% zScintA      Nt-by-1 vector containing the normalized complex scintillation
%              time history in the form of averages over Ts with sampling
%              interval Ts.  zScintA(kp1) is the average over tk to tkp1.
% 
% tScint       Nt-by-1 vector of time points corresponding to zScint
%
%+------------------------------------------------------------------+
% References: Humphreys et al., "Simulating Ionosphere-Induced Scintillation
% for Testing GPS Receiver Phase Tracking Loops," In preparation for
% submission to an IEEE journal.
% 
% Notes:  K is related to S4 under the Ricean model by 
%         K = sqrt(m^2 - m)/(m - sqrt(m^2 - m)) where 
%         m = max(1,1/(S4^2)).
%
% Author:  Todd Humphreys
%+==================================================================+


%----- Model parameters
TsSub = Ts/Nspa;
% For convenience, make Ns even
Ns = Nt*Nspa;
tScint = [0:Ns-1]'*TsSub;
df = 1/tScint(end);
beta0 =  1.23964643681047;

%----- Set up the Butterworth filter
% See Eq. 3.8 in the scintillation modeling paper
Bd = beta0/(sqrt(2)*pi*tau0);  
fn = 0.5*(1/TsSub);
Wn = Bd/fn;
[B,A] = butter(2,Wn);

%----- Filter white noise to produce a time history with the correct
%      power spectrum.
nVec = randn(Ns,1) + j*randn(Ns,1);
xiScint = filter(B,A,nVec);

%----- Add the line of sight component to get zScint
% Calculate the mean squared value sigmaxi2
sigmaxi2 = 0.5*mean(xiScint.*conj(xiScint));
% Calculate the required value for the line-of sight component zBar.
% This comes from the definition of the Ricean K parameter: K = (zBar^2/2)/sigmaxi2
zBar = sqrt(2*sigmaxi2*K);
% zScintC represents the high sampling rate (quasi continuous time) time
% history of the complex scintillation
zScintC = zBar + xiScint;
% Normalize zScintC so that E[|zScintC(k)|^2] = 1.
zScintC = zScintC/sqrt(mean(abs(zScintC).^2));

%----- Sample zScintC to get zScint
iidum = [1:Nspa:Ns]';
zScint = zScintC(iidum);
tScint = [0:Nt-1]'*Ts;

%----- Average zScintC to get zScintA
M = zeros(Nspa,Nt);
M(:) = zScintC;
zScintA = conj(mean(M,1)');

