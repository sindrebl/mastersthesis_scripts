% This script is used to process data extracted from ABAQUS using the
% funtion "writeOutputTofile" in "myUtilities.py". It calculates the
% acoustic impedance of the FE model and stores it in "impedanceVector". It
% also calculates the acoustic impedance using the Three Phase Method (TPM)
% and stores it in impedanceTPMVector. It is importent to notice that the
% materials used here are equal to those used in ABAQUS.
% Note that different simulations requires different methods for post
% processing. This is just a very general one used. 

clc;clear all;close all;
addpath('F:/masteroppgave/script/scripts-tobe-published/MATLAB/supportingFunctions/');

% script setup
script = {};
script.select_folder        = true;
script.resample             = true; 
script.resample_Ts          = false;

%========================================================================== 
%%  Choose directory and read infoFile
%==============================================================================
if script.select_folder
    fpath = uigetdir(pwd,'Select directory with ABAQUS results');
else
  fpath = 'F:/masteroppgave/script/textOutDirectory/Comparison/';
end

%==========================================================================
%% Define materials
%==========================================================================
silver.rho      = 10.49e-12;
silver.E        = 76.0e+6;
silver.nu       = 0.37;
PMMA.rho        = 1.16e-12;
PMMA.E          = 5.37e+06;
PMMA.nu         = 0.343;
epoxy.rho       = 1.12e-12;
epoxy.E         = 4.35e+06;
epoxy.nu        = 0.368;
polystyrene.rho = 1.05e-12;
polystyrene.E   = 3.45e+06;
polystyrene.nu  = 0.358;

%==========================================================================
%% Choose materials
%==========================================================================
coreMaterial = polystyrene;
shellMaterial = silver;
matrixMaterial = epoxy;

%==========================================================================
%% read infoFiles
%==========================================================================
[infoLines] = textread(fullfile(fpath, 'infoFile.txt'),'%s');
splitStr    = ',';
infoValues  = regexp(infoLines,splitStr,'split');

Nsims = numel(infoLines);

%==========================================================================    
%%  Read simulation files
%==========================================================================
sim = cell(1,Nsims);
for ss = 1:Nsims
    sim{ss}.infoLine = infoLines{ss};
    sim{ss}.name = char(infoValues{ss}(1));
    sim{ss}.NoOfElem = str2double(char(infoValues{ss}(5)));
    
    sim{ss}.radius = str2double(strrep(sim{ss}.name(strfind(sim{ss}.name,'R')+1:strfind(sim{ss}.name,'ST')-1),'_','.'));
    sim{ss}.shellT = str2double(strrep(sim{ss}.name(strfind(sim{ss}.name,'ST')+2:strfind(sim{ss}.name,'PC')-1),'_','.'));
    sim{ss}.PC = str2double(strrep(sim{ss}.name(strfind(sim{ss}.name,'PC')+2:strfind(sim{ss}.name,'N')-1),'_','.'));	

    sim{ss}.sT  = dlmread(fullfile(fpath,['tStr__',sim{ss}.name,'.txt']));
    sim{ss}.sB  = dlmread(fullfile(fpath,['bStr__',sim{ss}.name,'.txt']));
    sim{ss}.d   = dlmread(fullfile(fpath,['dist__',sim{ss}.name,'.txt']));
    sim{ss}.t   = dlmread(fullfile(fpath,['time__',sim{ss}.name,'.txt']));
    sim{ss}.w   = dlmread(fullfile(fpath,['width_',sim{ss}.name,'.txt']));
end

elemVector              = zeros(Nsims,1);
velocityVector          = zeros(Nsims,1);
impedanceVector         = zeros(Nsims,1);
impedanceTPMVector      = zeros(Nsims,1);
shellFracVector         = zeros(Nsims,1);

%==========================================================================
%% Study individual simulations
%==========================================================================
for ss = 1:Nsims
    % Rename variables
    t_s                 = sim{ss}.t;                % Time vector [s]
    [Nt,Nelem]          = size(sim{ss}.sT);         % Number of time sample, Number of elements
    simName             = sim{ss}.name;				% Simulation name
    Ts                  = t_s(2)-t_s(1);            % Sampling interval [s]
    
    sT                  = sim{ss}.sT;				% Stress in the top elements
    sB                  = sim{ss}.sB;				% Stress in the bottom elements
    d                   = sim{ss}.d;				% Distance between top and bottom
    
    radius 				= sim{ss}.radius;			% Core radius
    shellThickness 		= sim{ss}.shellT;			% Shell thickness
    partConc			= sim{ss}.PC;				% particle concentration
    
    volumeSphere = (4.0/3)*pi*(radius+shellThickness)^3;
    volumeCore = (4.0/3)*pi*radius^3;
    volumeShell = volumeSphere-volumeCore;
    width =  (volumeSphere/(2*partConc*pi))^(1 / 3.0); 
    cellHeight = width;
    volumeCell = pi*(width^2)*2*cellHeight;

    sphere.rho = density_eff( (volumeCore/volumeSphere),coreMaterial.rho,shellMaterial.rho);
    composite.rho = density_eff(partConc,sphere.rho,matrixMaterial.rho);        
    shellFraction =(1-(volumeCore/volumeSphere))*partConc;
    sim{ss}.shellFrac = shellFraction;
    sim{ss}.composite_rho = composite.rho;

%==========================================================================   
%%  Resample
%==========================================================================
    if script.resample
        %resample using number of points
        sim{ss}.t_rs_s = linspace(0,t_s(end),Nt)';          % new time vector 
        Ts = sim{ss}.t_rs_s(2)-sim{ss}.t_rs_s(1);           % Sampling interval [s]
        
        sim{ss}.sB_rs = interp1(t_s,sim{ss}.sB,sim{ss}.t_rs_s,'spline',0);
        sim{ss}.sT_rs = interp1(t_s,sim{ss}.sT,sim{ss}.t_rs_s,'spline',0);
         
        %% Rename variables
        t_s_old             = t_s;
        t_s             	= sim{ss}.t_rs_s;
        [Nt,Nelem]          = size(sim{ss}.sT_rs);
        sT       			= sim{ss}.sT_rs;
        sB                  = sim{ss}.sB_rs;
    end
    
%==============================================================================    
%%  prepend/append zeros - used to do fourier analysis during development
%============================================================================== 
    Nzero_app = 1000;
    Nzero_pre = 1000;
    
    sT_padd = [zeros(Nzero_pre,Nelem); sT; zeros(Nzero_app,Nelem)];
    sB_padd = [zeros(Nzero_pre,Nelem); sB; zeros(Nzero_app,Nelem)];
    
    Nt_padd = size(sT_padd,1);
    t_padd  = [(-Nzero_pre-1:-1)*Ts, t_s', t_s(end)+(1:Nzero_app-1)*Ts]';

    f_plot_padd = (0:Nt_padd-1)/(Nt_padd*Ts); 				% Frequency vector for plotting [Hz]
    f_calc_padd = (-Nt_padd/2:Nt_padd/2-1)'/(Nt_padd*Ts);  	% Frequency vector for calculations [Hz]
    
%==============================================================================
%% Apply window on measurement
%==============================================================================
    windowp.pulselength_s = 3.5*1e-7;  % 55.4-54.4 µs
    windowp.tukey_slopeprcent = 0.30;
    windowp.windowlength = ceil((windowp.pulselength_s/Ts) / (1 - windowp.tukey_slopeprcent));

    sT_padd_win = zeros(size(sT_padd));
    sB_padd_win = zeros(size(sT_padd));

    sT_padd_old = sT_padd;
    sB_padd_old = sB_padd;
    
    for el_ii = 1:Nelem
        % Window on Top
        [~, sT_maxidx] = max(sT_padd(:,el_ii));
        win_start_idx   = sT_maxidx - floor(windowp.windowlength/2)+1;
        win_end_idx     = sT_maxidx + ceil(windowp.windowlength/2);
        
        totalwindow_T = [ zeros(win_start_idx-1,1);...
                            tukeywin(windowp.windowlength,windowp.tukey_slopeprcent);...
                            zeros(Nt_padd-win_end_idx,1)];
        sT_padd_win(:,el_ii) = sT_padd(:,el_ii).*totalwindow_T;


        % Window on Bottom
        [~, sB_maxidx] = max(sB_padd(:,el_ii));
        win_start_idx   = sB_maxidx - floor(windowp.windowlength/2)+1;
        win_end_idx     = sB_maxidx + ceil(windowp.windowlength/2);

        totalwindow_B = [ zeros(win_start_idx-1,1);...
                            tukeywin(windowp.windowlength,windowp.tukey_slopeprcent);...
                            zeros(Nt_padd-win_end_idx,1)];
        sB_padd_win(:,el_ii) = sB_padd(:,el_ii).*totalwindow_B;      
    end    
%==============================================================================
%% XCORR analysis
%==============================================================================
    for elem = 1:Nelem
        sim{ss}.tau_xcorr_win(elem) = Ts*xcorr_shift(sB_padd_win(:,elem),sT_padd_win(:,elem));
    end
    sim{ss}.velocity = d./sim{ss}.tau_xcorr_win';
    
%============================================================================
%% weighted velocity
%============================================================================
    elemVector(ss) = sim{ss}.NoOfElem;
    totWidth = sum(sim{ss}.w);
    velocityWeighted  = 0;
    radius = 0;
    for elem = 1:length(sim{ss}.w)
        velocityWeighted = (velocityWeighted ...
        + sim{ss}.velocity(elem).*((radius+sim{ss}.w(elem))^2 ...
        - (radius^2))/ totWidth^2 );
        radius = radius + sim{ss}.w(elem);
    end
    velocityVector(ss) = velocityWeighted*1e-6;  %[m/s]
    impedanceVector(ss) = velocityWeighted*sim{ss}.composite_rho*1e3; %[Mrayl]
    shellFracVector(ss) = sim{ss}.shellFrac;
    compoTPM = homogenization_threePhase(coreMaterial,shellMaterial,matrixMaterial,(sim{ss}.radius^3/((sim{ss}.radius+sim{ss}.shellT)^3)),sim{ss}.PC);
    impedanceTPMVector(ss) = compoTPM.imp*1000; % [Mrayl]
end

disp(impedanceVector);
disp(velocityVector);



