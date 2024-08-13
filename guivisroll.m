function [] = guivis()
% GUIVISROLL   Graphical User Interface for visualizing scintillation,
%              rolling version
%
%+------------------------------------------------------------------+
% References:
%
% Notes: The scintillation data zkhist and tkhist are assumed to be found
%        in the file scintDat.mat.
%
% Author:  Todd Humphreys, Aug 2, 2007
%+==================================================================+


%----- Parameters
tWin = 20;        % Window length (sec)
dtRoll = 0.15;       % Roll increment
pauseTime = 0.01;
% Detrending filter cutoff
fD = 0.2;

%----- Local variables
S = load('scintDat');
akhist = abs(S.zkhist);                   % amplitude
pkhist = unwrap(angle(S.zkhist))/2/pi;    % phase in cycles
tkhist = S.tkhist;
mScint = mean(log(akhist));
Ts = tkhist(2) - tkhist(1);
dii = max(1,floor(dtRoll/Ts));
Nii = length(tkhist);
NWii = floor(tWin/Ts);
ii = 1;
iiMax = Nii - NWii;
if(iiMax<0)
    warning('Data length must be at least 20 secs. long to plot, exiting.')
    return
end
% rollFlag = -1 (roll left), 0 (stop), 1 (roll right).
rollFlag = 1; 

%----- High-pass filter the phase data
fs = 1/Ts;        % sampling frequency (Hz)
Wn = fD/(0.5*fs);
[B,A] = butter(4,Wn,'low');
pkhistL = filtfilt(B,A,pkhist);
pkhist = pkhist - pkhistL;

f = figure('visible','off','name', 'Scintillation Visualizer');
hAxis = axes('position', [0.1300    0.2100    0.7750    0.7150]);

hPrev  = uicontrol('Style','pushbutton',...
                  'String','<-Roll','Position',[45,5,60,25],...
                  'callback', {@prevCallback});

hNext  = uicontrol('Style','pushbutton',...
                  'String','Stop','Position',[480,5,60,25],...
                  'callback', {@nextCallback});

hClose = uicontrol('Style','pushbutton',...
                  'String','Close','Position',[265,5,60,25],...
                  'callback', {@closeCallback});


%----- Normalize units
handleVec = [hAxis,hPrev,hNext,hClose];
set(handleVec,'units','normalized');

set(f,'visible','on');

%----- Initial plot
iidum = [1:NWii]';
cla;
hP1 = plot(tkhist(iidum),log(akhist(iidum)), 'linewidth', 2);
hold on;
hP2 = plot(tkhist(iidum),(pkhist(iidum)) + mScint, 'g', 'linewidth', 2);
grid on; ylim(round(mScint) + [-4 2]); xlabel('Time (s)');
ylabel('LOG(\alpha(t_k)) and \theta(t_k) (cycles)');
legend('Amplitude', 'Phase');

rollScint

%----- Callback functions
function prevCallback(hObject,eventdata)
str = get(hObject,'string');
if(strcmp(str,'<-Roll'))
  set(hObject,'string','Stop');
  set(hNext,'string','Roll->');
  rollFlag = -1;
else
  set(hObject,'string','<-Roll');
  set(hNext,'string','Roll->');
  rollFlag = 0;
end
rollScint
end

function nextCallback(hObject,eventdata)
str = get(hObject,'string');
if(strcmp(str,'Roll->'))
  set(hObject,'string','Stop');
  set(hPrev,'string','<-Roll');
  rollFlag = 1;
else
  set(hObject,'string','Roll->');
  set(hPrev,'string','<-Roll');
  rollFlag = 0;
end
rollScint
end

function closeCallback(hObject,eventdata)
close(f);
rollFlag = 0;
end

%----- The rollScint function
function rollScint
while(rollFlag)
  ii = ii + rollFlag*dii;
  ii = max(1,ii);  ii = min(iiMax,ii);
  iidum = [ii:1:NWii + ii-1]';
  if(ishandle(hP1)&ishandle(hP2))
  set(hP1,'ydata', log(akhist(iidum)));
  set(hP2,'ydata', (pkhist(iidum)) + mScint);
  drawnow
  pause(pauseTime);
  else
    return;
  end
end
end

end
