function guiscint
% GUISCINT   Graphical User Interface for the scintillation simulator
%
% INPUTS
% input1       Definition
%
% OUTPUTS
% 
%+------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys, 
%+==================================================================+

close all;
%----- Local parameters
% Nspb is the number of sub-samples per accumulation
Nspa = 8; % Number of sub-samples that are used in calculating the averages in zScintA  
fnh = 'helvetica';
fnt = 'times new roman';
Tb = 0.02;
% Initial values
ivS4 = 0.5; ivtau0 = 1; ivC_N0 = 45; ivTsim = 300;
% Extreme values
minS4 = 0; maxS4 = 1;
mintau0 = -1; maxtau0 = 0.3010;  % log10 of tau0
minC_N0 = 35; maxC_N0 = 55;      % dB-Hz
minTsim = 10; maxTsim = 1000;     % seconds
% The range of Te that will be plotted (logarithmically) on the graph
TeLow = 0.626; TeHigh = 3600*20; 
PeLow = Tb/TeHigh; PeHigh = Tb/TeLow;
PeLowLog = log(PeLow); PeHighLog = log(PeHigh);
Levels = {'Negligible','Weak', 'Moderate', 'Severe', ...
          'Catastrophic'};
NLevels = length(Levels);


%----- Initialize and hide the GUI as it is being constructed.
f =figure('Visible','off','menubar','none','Position',[360,500,396,320]);


%----- Construct the console frames
xConsole = 230; wConsole = 154; xTitles = 233;
hConsole1 =uicontrol( ...
    'Style','frame', ...
    'Position',[xConsole,170,wConsole,146], ...
    'BackgroundColor',[0.50 0.50 0.50]);
htConsole1  = uicontrol('style','text','string','Scintillation Parameters',...
                        'Position',[xTitles,300,145,10], ...
                        'BackgroundColor',[0.50 0.50 0.50],... 
                        'foregroundColor',[0,0,0.5],'fontweight','bold',...
                        'fontname',fnh);
hConsole2 =uicontrol( ...
    'Style','frame', ...
    'Position',[xConsole,102,wConsole,64], ...
    'BackgroundColor',0.5*[1,1,1]);
htConsole2   = uicontrol('style','text','string','Expected Nominal C/N0',...
                         'Position',[xTitles,152,145,12], ...
                         'BackgroundColor',0.5*[1,1,1],... 
                         'foregroundColor',[0,0,0.5],'fontweight','bold',...
                         'fontname',fnh);

hConsole3 =uicontrol( ...
    'Style','frame', ...
    'Position',[xConsole,4,wConsole,94], ...
    'BackgroundColor',0.5*[1,1,1]);
htConsole3   = uicontrol('style','text','string','Simulation Settings',...
                         'Position',[xTitles,84,145,12], ...
                         'BackgroundColor',0.5*[1,1,1],... 
                         'foregroundColor',[0,0,0.5],'fontweight','bold',...
                         'fontname',fnh);

%----- Construct the components.
xSlider = 244;wSlider = 125;hSlider = 15;wEdit = 50;xEdit = 279;
hS4      = uicontrol('Style','slider','min',minS4,'max',maxS4,...
                     'value',ivS4,...
                     'position',[xSlider,260,wSlider,hSlider],...
                     'sliderstep',[0.05,0.1],...
                     'callback', {@S4sliderCallback});
ivtau0Log10 = log10(ivtau0);
htau0    = uicontrol('Style','slider','min',mintau0,'max',maxtau0,...
                     'value',ivtau0Log10,...
                     'position',[xSlider,200,wSlider,hSlider],...
                     'sliderstep',[0.05,0.1],...
                     'callback', {@tau0sliderCallback});
hC_N0    = uicontrol('Style','slider','min',minC_N0,'max',maxC_N0,...
                     'value',ivC_N0,...
                     'position',[xSlider,130,wSlider,hSlider],...
                     'sliderstep',[0.05,0.1],...
                     'callback', {@C_N0sliderCallback});
htS4     = uicontrol('style','text','string','S4 Index',...
                     'Position',[330-56,280,60,10], ...
                     'BackgroundColor',[0.50 0.50 0.50],... 
                     'foregroundColor',[0,0,0],'fontweight','bold',...
                     'fontname',fnh);
httau0   = uicontrol('style','text','string','tau0',...
                     'Position',[274,220,60,10], ...
                     'BackgroundColor',[0.50 0.50 0.50],... 
                     'foregroundColor',[0,0,0],'fontweight','bold',...
                     'fontname',fnh);
httau0U  = uicontrol('style','text','string','sec',...
                     'Position',[332,183,20,12], ...
                     'BackgroundColor',[0.50 0.50 0.50],... 
                     'foregroundColor',[0,0,0],'fontname',fnh);
heS4     = uicontrol('style','edit','position',[xEdit,243,wEdit,hSlider],...
                     'value',ivS4,'string',num2str(ivS4),...
                     'callback', {@S4editCallback});
hetau0   = uicontrol('style','edit','position',[xEdit,183,wEdit,hSlider],...
                     'value',ivtau0,'string',num2str(ivtau0),...
                     'callback', {@tau0editCallback});
heC_N0   = uicontrol('style','edit','position',[xEdit,113,wEdit,hSlider],...
                     'value',ivC_N0,'string',num2str(ivC_N0),...
                     'callback', {@C_N0editCallback});
htC_N0U  = uicontrol('style','text','string','dB-Hz',...
                     'Position',[332,113,35,12], ...
                     'BackgroundColor',[0.50 0.50 0.50],... 
                     'foregroundColor',[0,0,0],'fontname',fnh);

hTa = uicontrol('Style','popupmenu',...
           'String',{'10','20'},...
           'Position',[244,60,40,20], 'value', 1);
htTa   = uicontrol('style','text','string','Ta (ms)',...
                   'Position',[284,65,50,12], ...
                   'BackgroundColor',[0.50 0.50 0.50],... 
                   'foregroundColor',[0,0,0],'fontweight','bold',...
                   'fontname',fnh);
heTsim   = uicontrol('style','edit','position',[244,25,40,hSlider],...
                   'value',ivTsim,'string',ivTsim,... 
                   'callback', {@TsimeditCallback});
hteTsim   = uicontrol('style','text','string','Length (sec)',...
                   'Position',[290,25,70,12], ...
                   'BackgroundColor',[0.50 0.50 0.50],... 
                   'foregroundColor',[0,0,0],'fontweight','bold',...
                   'fontname',fnh);
hSim = uicontrol('Style','pushbutton',...
                 'String','Simulate','Position',[45,5,60,25],...
                 'callback', {@SimCallback});
htSim   = uicontrol('style','text','string',...
                    'Data written to scintDat.mat',...
                    'Position',[120,4,80,25], ...
                    'BackgroundColor',[0.8 0.8 0.8],... 
                    'foregroundColor',[0,0,0],...
                    'fontname',fnh, 'visible', 'off');

%----- The axis and its labels
hAxis = axes('Units','pixels','Position',[50,40,50,260],...
             'visible', 'off');
htAxis = uicontrol('style','text','string','Te = 10 sec',...
                   'Position',[35,300,80,12], ...
                   'foregroundColor',[0,0,0.5],...
                   'backgroundColor',0.8*[1,1,1],...
                   'fontname',fnh,'visible', 'off');
hLevelVec = zeros(NLevels,1);
ybase = 33;
dy = 62;
for ii=1:NLevels
  str = Levels{ii};
  ylab = ybase + dy*(ii-1);
  M = hot;
  clr = M((ii-1)*6 + 1,:);
  hLevelVec(ii) = uicontrol('style','text','string',str,...
                            'Position',[103,ylab,80,20], ...
                            'foregroundColor',clr,...
                            'backgroundColor',0.8*[1,1,1],...
                            'fontname',fnh,...
                            'fontweight','bold','fontsize',10,...
                            'horizontalalignment', 'left');
end

%align([htConsole1,httau0],'Center','None');

%----- Initialize the GUI.
% Change units to normalized so components resize automatically.
handleVec = [f,hConsole1,htConsole1,hConsole2,htConsole2,hConsole3, ...
             htConsole3,hS4,htau0,hC_N0,htS4,httau0,heS4,hetau0,heC_N0, ...
             hTa,htTa,heTsim,hteTsim,hSim,hAxis,httau0U,htC_N0U,...
             htAxis,htSim,hLevelVec'];
set(handleVec,'Units','normalized');
% Assign the GUI a name to appear in the window title.
set(f,'Name','Equatorial Scintillation Simulator')
movegui(f,'northeast')
% Make the GUI visible.
set(f,'Visible','on')
updateGraph

%----- Program the Callbacks
function S4sliderCallback(hObject,eventdata) 
v = get(hObject,'value');
v = fix(100*v)/100;
set(heS4,'value',v,'string',num2str(v));
updateGraph
end

function tau0sliderCallback(hObject,eventdata) 
v = get(hObject,'value');
v = 10^v; v = round(100*v)/100;
set(hetau0,'value',v,'string',num2str(v));
updateGraph
end

function C_N0sliderCallback(hObject,eventdata)
v = get(hObject,'value');
set(heC_N0,'value',v,'string',num2str(v));
updateGraph
end

function S4editCallback(hObject,eventdata)
vstring = get(hObject,'string');
v = eval(vstring);
v = fix(100*v)/100;
if(v < minS4) v = minS4; end
if(v > maxS4) v = maxS4; end
set(hS4,'value',v);
set(heS4,'value',v, 'string', num2str(v));
updateGraph
end

function tau0editCallback(hObject,eventdata)
vstring = get(hObject,'string');
v = (eval(vstring));
v = round(100*v)/100; v = log10(v);
if(v < mintau0) v = mintau0; end
if(v > maxtau0) v = maxtau0; end
v = 10^v;
set(hetau0,'value',v,'string',num2str(v));
set(htau0,'value',log10(v));
updateGraph
end

function C_N0editCallback(hObject,eventdata) 
vstring = get(hObject,'string');
v = eval(vstring);
v = round(10*v)/10;
if(v < minC_N0) v = minC_N0; end
if(v > maxC_N0) v = maxC_N0; end
set(hC_N0,'value',v);
set(heC_N0,'value',v,'string',num2str(v));
updateGraph
end


function TsimeditCallback(hObject,eventdata) 
vstring = get(hObject,'string');
v = fix(eval(vstring));
v = round(v);
if(v < minTsim) v = minTsim; end
if(v > maxTsim) v = maxTsim; end
set(heTsim,'string',num2str(v),'value',v);
end

function SimCallback(hObject,eventdata) 
set(htSim,'visible','off');
str = get(hSim,'string');
if(strcmp(str,'Simulate'))
  set(hSim,'string', 'Busy', 'enable', 'off');
  pause(0.01);
  v = get(hTa,'value'); vstr = get(hTa,'string');
  Ta = eval(vstr{v})/1000; % Sampling interval of the output time histories zScint and zScintA, in seconds.
  tau0 = get(hetau0,'value'); % decorrelation time of the second-order Butterworth model (first input)
  S4 = get(heS4,'value'); m = max(1,1/(S4^2)); % scintillation index (second input)
  K = sqrt(m^2 - m)/(m - sqrt(m^2 - m)); % Eq. (5)
  Tsim = get(heTsim,'value'); % simulation time
  Napb = Tb/Ta; % ratio of Sampling interval of the output time histories zScint and zScintA and the bit sampling time
  Nt = Napb*round(Tsim/Tb); % Number of samples in the output time histories zScint and zScintA. = % of samples used in time histories zScint and zScintA * total samples
  [zskhist,zkhist,tkhist] = scintModel04(Ta,Nt,tau0,K,Nspa);
  save scintDat zkhist tkhist;  
  set(htSim,'visible','on');
  set(hSim,'string', 'View', 'enable', 'on');
else
  set(hSim,'string', 'Simulate', 'enable', 'on');
  guivisroll   % or guivis
end

end

%----- The update function
function updateGraph
set(hSim,'string', 'Simulate', 'enable', 'on');
set(htSim,'visible','off');
S4 = get(heS4,'value');
tau0 = get(hetau0,'value');
C_N0 = get(heC_N0,'value');
tau0Vec = tau0*ones(3,1);
PeVec = estimatePe(S4,tau0Vec,C_N0);
PeLog = log(PeVec(2) + eps);
if(PeLog > PeHighLog)
  error('Scintillation severity out of range');
end
if(PeLog < PeLowLog) PeLog = PeLowLog; end
yVal = (PeLog - PeLowLog)/(PeHighLog-PeLowLog);
hbar = bar(yVal);
set(hAxis,'ylim',[0 1],'xtick',[],'ytick',[],...
          'xticklabel',[],'yticklabel',[], 'visible', 'on');
Te = Tb/exp(PeLog);
if(Te <= 10)
  Te = round(10*Te)/10;
  ostr = ['Te = ' num2str(Te) ' sec'];
elseif(Te <= 1000)
  Te = round(Te);
  ostr = ['Te = ' num2str(Te) ' sec'];
elseif((log(PeVec(2)+eps) < PeLowLog))
  Te = TeHigh;
  Te = round(10*Te/3600)/10;
  ostr = ['Te > ' num2str(Te) ' hrs'];
else  
  Te = round(10*Te/3600)/10;
  ostr = ['Te = ' num2str(Te) ' hrs'];
end
set(htAxis,'string',ostr, ...
           'visible', 'on');
end

end