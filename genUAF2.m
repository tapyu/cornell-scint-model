function [] = genUAF2(zkhist,tkhist,prn,Tquiet)
% genUAF2        Generate a User Actions File from a scintillation time history
%
% Created by Joanna Hinks and Todd Humphreys, 29 May 2008.
%
% INPUTS
% zkhist         Nt-by-2-by-Ns matrix containing the normalized complex scintillation
%                time histories that will be used to drive variations in the
%                L1 and/or L2 output signal of the RF signal simulator.
%                Column 1 contains the L1 scintillation history; column 2
%                contains the L2 scintillation history. A column of zeros
%                means no scintillation will be commanded for that
%                frequency. The dimension Ns corresponds to the satellite 
%                prn number. The time history is expressed in the form of 
%                averages over Ts with sampling interval Ts.  zkhist(kp1,L,s) 
%                is the average over tk to tkp1 for the frequency L and prn s.
%
% tkhist         Nt-by-1 vector of time points corresponding to zkhist.
%                The sampling interval Ts = tkhist(ii+1)  - tkhist(ii). 
%
% prn            Ns-by-1 vector of PRN identifiers of the GPS satellites 
%                whose signals will be commanded to scintillate (in the same 
%                order in which the 3-dimensional zkhist matrix is
%                stacked).
%
% Tquiet         Length of the quiet interval (no scintillation) that
%                precedes the onset of the scintillation interval, in
%                seconds.  The quiet interval allows the tracking loops of 
%                the receiver under test to lock and settle.
%

% Check dimensions of zkhist (must be Nt x 2 x Ns)
if length(zkhist(1,:,1)) ~= 2
    error('Input zkhist is not correctly formatted.')
end
Ns = length(prn);
if length(zkhist(1,1,:)) ~= Ns
    error('Mismatch of signals and PRNs.');
end
    
% Prompt user for output file name.  Use the file extension *.cmd, although
% the resulting files can be opened as text files.
output_file = input('Please enter output file name: ','s');

% Determine constant sampling rate (in units of csec)
Ts = 100*(tkhist(2)-tkhist(1));
Nt = length(tkhist);

% Determine Tonset and parse it into units of day,hour,min,sec,csec
Tonset = Tquiet+Ts/100;
Tonset = round(Tonset*100)/100;
day = floor(Tonset/(60*60*24)); Tresidual = Tonset-day*60*60*24;
hour = floor(Tresidual/(60*60)); Tresidual = Tresidual-hour*60*60;
min = floor(Tresidual/60); Tresidual = Tresidual-min*60;
sec = floor(Tresidual); Tresidual = Tresidual-sec;
csec = round(Tresidual*100);

% Determine which frequencies will be commanded to scintillate
noL1 = 0; noL2 = 0;
if sum(zkhist(:,1,:)~=0) == 0
    noL1 = 1;
end
if sum(zkhist(:,2,:)~=0) == 0
    noL2 = 1;
end

if noL1 && noL2
    error('No scintillation commanded on any frequency');
elseif noL1
    zkhist = zkhist(:,2,:);
elseif noL2
    zkhist = zkhist(:,1,:);
end

% Set up parameters to keep track of which frequencies are present
NL = length(zkhist(1,:,1));
if NL == 2
    Lvec = [0 1];
else
    Lvec = noL1; 
end


% Separate out amplitude and phase information. 
amphist = abs(zkhist);
phasehist = unwrap(angle(zkhist));

% Compute change in signal level from amphist
sig_level = 10*log10(amphist.^2);

% Convert change in phase from radians to meters
lambdavec = [0.190293672798365,0.244210213424568];
for ii = 1:NL
    phasehist(:,ii,:) = (phasehist(:,ii,:)/(2*pi))*lambdavec(Lvec(ii)+1);
end


% Create file and write command lines
fid = fopen(output_file,'w');
fprintf(fid,'NBLK2\n');
for k = 1:Nt
    for f = 1:NL
        for s = 1:Ns
            output = sprintf('%1d %02d:%02d:%02d.%02d,MOD,v1_a1,gps,%02d,0,0,0,%1d,0,%2.3f,%2.3f,0\n',day,hour,min,sec,round(csec),prn(s),Lvec(f),sig_level(k,f,s),phasehist(k,f,s));  
            fprintf(fid,output);
        end
    end
    
    % Update timestamp to advance by Ts
    csec = csec + Ts;
    if (csec >= 100)
        sec = sec + 1;
        csec = csec - 100;
    end
    if(sec >= 60)
        min = min + 1;
        sec = 0;
    end
    if(min >= 60)
        hour = hour+1;
        min = 0;
    end;
    if(hour >= 24)
        day = day+1;
        hour = 0;
    end;
end

% Zero the output once scintillation has ended
for f = 1:NL
    for s = 1:Ns
        output = sprintf('%1d %02d:%02d:%02d.%02d,MOD,v1_a1,gps,%02d,0,0,0,%1d,0,0,0,0\n',day,hour,min,sec,csec,prn(s),Lvec(f));
        fprintf(fid,output);
    end
end

% Close file
fclose(fid);

