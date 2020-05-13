%% Formin tracks pipeline
%
% This program is used to measure the speed of formins
% It needs as an input a .tif file extracted from a movie and the creation 
% of an Image Sequence in a fov1 folder, these images (files) must be named 
% fovn_ (n = embryo number) 
% create a folder containing all the fovn folders with the image sequences
%
%
% Step 1 : define the number of embryos (n), tracking and Plot D Alpha 
% parameters. 
% Step 2 : Create a Kilfoil Work Directory. Create a unique time matrix in 
% the kilfoil folder and copy/past every fovn folder as a fov1 folder 
% containing fov1_ files so it can be the input of Kilfoil_Pretrack_GUI_v2, 
% run Kilfoil_Pretrack_GUI_v2 and save the chosen values of MassCut and 
% Imin before closing it.
% Step 3 : Tracking routine : mpretrack and fancytrack, create and save 
% the matrix resEmb with the parameters of every track.
% Step 4 : Fill gaps and Dedrift.
% Step 5 : Choose the values Max and Min of the Sup- and Sub- populations 
% ans display Plot D Alpha.
% Step 6 : Display the plots of every track and pick the plots with 
% linear tracks then measure the lifetime and the length of every track
% and then the speed of every formin and the mean speed of the formins 
% by embryo.
%
%
%
% Last modified 10-May-2020 19:56 Anais DJEBARA
%% Define the number of embryos and enter tracking parameters

prompt = 'Number of embryos = ' ;
num_of_embryos = input(prompt);

prompt = 'MovieTimeStep value = ' ;
MovieTimeStep = input(prompt);

featsize = 3 ; barrRg = 50 ; barrCc = 1 ; IdivRg = 0 ; field = 2 ; 
maxdisp = 2 ; goodenough = 1 ; memory = 2 ; barrI = 1 ; 

alphaSup = 1.1 ; alphaSub = 1.1 ; % Plot D Alpha


 for i = 1:num_of_embryos   
%% Create a Kilfoil Work Directory and copy the fovn files

     mkdir('KilfoilWorkDir'); % Kilfoil Work Directory
     
     copyfile(['fov', num2str(i)], 'KilfoilWorkDir/fov1'); % copy the fovn file
     cd('./KilfoilWorkDir/');   
    
     time=[0:.5:10000];
     save(['fov1_times.mat'],'time');
     % creates and saves a time matrix for the Kilfoil_Pretrack_GUI
   
    if i ~= 1 
        cd('./fov1')
        content = dir('*fov*');
        NbFiles = length(content)-1;
        for j = 0:NbFiles
            copyfile(['fov', num2str(i), '_', sprintf('%04d',j), '.tif'], ['fov1_', sprintf('%04d',j), '.tif']);
            delete(['fov', num2str(i),'_', sprintf('%04d',j), '.tif']);
        end
    % copies the fovn_ files and pasts them as fov1_ files in the KilfoilWorkDir file
    end
    Imin = 0 ; MassCut = 0 ; 
     cd('../')
     Kilfoil_Pretrack_GUI_v2
     % a version of KilFoil_Pretrack_GUI that allows to save an Imin and 
     % a MassCut value in the workspace
     
     fig = uifigure;
     fig.Position = [500 500 500 350];
     uialert(fig,'Press Ok or Close me only if you are done and have pressed saved in Kilfoil GUI.', ...
    'Program Information','Icon','info','CloseFcn','uiresume(fig)')
     uiwait(fig)
     close(gcf); 
     % save the values of MassCut and Imin before closing Kilfoil_Pretrack_GUI_v2
     
%%  Tracking routine : mpretrack / fancytrack / embEnv / resEmbn
     
    if i ~= 1
       cd('../')
    end 
    time=[0:.5:10000]; 
    save(['fov', num2str(i), '_times.mat'],'time'); 
    % saves a time matrix for every fovn used in mpretrack
    
    if Imin == 0 && MassCut == 0
        warning('Save MassCut and Imin values to continue')
    else
        
    nbframes = length(dir(['fov', num2str(i)]))-4 ; % frame number
        
    mpretrack('./',i,featsize,barrI,barrRg,barrCc,IdivRg,nbframes,MassCut,Imin,field) ;
    
    fancytrack('./',i, featsize, maxdisp, goodenough, memory ) ;
      
    end
    
    addpath(genpath('Bead_trackig'), 'Feature_finding') ; 
    % adds Bead-tracking and Feature_finding to matlab path
    cd('./Bead_tracking/res_files') ;
    load(['res_run',num2str(i),'.mat']); % load res_run(i).mat
    cd(['../../fov', num2str(i)]) ;
    
    embEnv = emb_env(['fov',num2str(i), '_0001.tif']); %%
    [embEnv(:,1),embEnv(:,2)] = poly2cw(embEnv(:,1),embEnv(:,2)); % draw a polygone around the embryo
    resEmb = res(inpolygon(res(:,1),res(:,2),embEnv(:,1),embEnv(:,2))==1,:);
    
    cd(['../'])
    outfname = strcat(['resEmb', num2str(i),'.mat']);
    save(outfname, 'resEmb');
    save(['resEmb', num2str(i),'.mat'], 'res') ;
    % create and save a res matrix
    
%%    Fill gaps and Dedrift

    kf = kf2sp_addgaps(res) ; % convert from Kilfoil to simple format (fill gaps with NaNs)
    kf_ = interpolateGaps(kf) ; % interpolate gaps
    res = sp2kf(kf_,MovieTimeStep) ; % convert from simple format to Kilfoil
    
    % dedrift 
    disp('Currently dedrifting')
    [resd,~] = dedrift(res, MovieTimeStep);
    % dedrifting of the input data to remove the effects of global drift 
    % from the mobility analysis, removes the effects of network's/embryo's 
    % drift from the mobility of the formins
 
    save(['resd', num2str(i) ,'.mat'],'resd'); % saves the new dedrifted res

%%  Choose the values Max and Min of the Sup- and Sub- populations and display Plot,D,Alpha
    
    % As we need to measure the speed of the moving formins only (superdiffusives) 
    % within two populations (Sub and Sup)
    
    load(['resd', num2str(i) ,'.mat'],'resd'); % load the dedrifted res
    kf = kf2sp_addgaps(resd); % fill gaps
    
    % define the parameters
    prompt = 'maxTauSup value = ' ; 
    maxTauSup = input(prompt);
    prompt = 'minTLSup value = ' ; 
    minTLSup = input(prompt);
    % max Tau used to compute the superdiffusive subpopulation
    % min Track Length used to compute the superdiffusive subpopulation
    
    prompt = 'maxTauSub value = ' ; 
    maxTauSub = input(prompt);
    prompt = 'minTLSub value = ' ; 
    minTLSub = input(prompt) ;
    % max Tau used to compute the subdiffusive subpopulation
    % min Track Length used to compute the subdiffusive subpopulation
       
    % plot D alpha function with tracks longer than minTLSup and interpolation window to maxTauSup
    [Alpha,D,ids] = plotD_Alpha(kf,maxTauSup,minTLSup,MovieTimeStep);
    
    superdif = ids(Alpha>alphaSup); % find the superdiffusive tracks (Alpha > alphaSup)
    subLongdif = ids(Alpha<alphaSup);
    
    % plot D alpha function with tracks longer than minTLSub and interpolation window to maxTauSub
    [Alpha2,D2,ids2] = plotD_Alpha(kf,maxTauSub,minTLSub,MovieTimeStep);
    
    superdif2 = ids(Alpha>alphaSub); % find the superdiffusive tracks (Alpha > alphaSub)
    subdif = unique(resd(~ismember(resd(:,8),superdif2),8));
    
    resSub = res(ismember(res(:,8),subdif),:);
    resSup = res(ismember(res(:,8),superdif),:); % populations of superdiffusive
    resSubLong = res(ismember(res(:,8),subLongdif),:); % corresponding population of subdiffusive
    
    prompt = 'Do you want to change your values? 1 if yes, 0 if no : ' ; 
    answer = input(prompt) ;
    if answer == 1
        while answer == 1 
            prompt = 'maxTauSup value = ' ; 
            maxTauSup = input(prompt);
            prompt = 'minTLSup value = ' ; 
            minTLSup = input(prompt);
            prompt = 'maxTauSub = ' ; 
            maxTauSub = input(prompt);
            prompt = 'minTLSub = ' ; 
            minTLSub = input(prompt) ;
            [Alpha,D,ids] = plotD_Alpha(kf,maxTauSup,minTLSup,MovieTimeStep);
            superdif = ids(Alpha>alphaSup); subLongdif = ids(Alpha<alphaSup);
            [Alpha2,D2,ids2] = plotD_Alpha(kf,maxTauSub,minTLSub,MovieTimeStep);
            superdif2 = ids(Alpha>alphaSub); 
            subdif = unique(res(~ismember(res(:,8),superdif2),8));
            resSub = res(ismember(res(:,8),subdif),:); 
            resSup = res(ismember(res(:,8),superdif),:);
            resSubLong = res(ismember(res(:,8),subLongdif),:);
            prompt = 'Do you want to change your values? 1 if yes, 0 if no : ' ; answer = input(prompt);
        end
    elseif answer == 0
    end

%%    Choose the plots with linear tracks / track length and formins' speed measurements
   
    kf = kf2sp_addgaps(res);
    [ Alpha,D, ids ] = plotD_Alpha(kf,13,20,MovieTimeStep);
    supids = ids(Alpha>1.5); % find the population with alpha > 1.5 
    % = superdiffusible population
    l = length(supids) ; 
    disp('Number of plots =') 
    disp(l)
    disp('press on space bar to keep the track , press on Enter key to exclude the track')
    count = 0;
    supidsLin = NaN;
   for j = supids
       figure;
       count = count+1;
       resSupLoc = [];
       resSupLoc = res(res(:,8)==j,:);
       plot(resSupLoc(:,1),resSupLoc(:,2));
       daspect([1 1 1]) % 1 c y et z
       testforsave = waitforbuttonpress();
       if testforsave == 0
       elseif testforsave == 1
           supidsLin = [supidsLin, j]; % saves the values of the picked plots = tracks
       end
       close(gcf) % closes and reopens a new plot
       % if mod(count,100) == 0,waitforbuttonpress();figure; count = 0;,end
   end
   supidsLin = supidsLin(2:end);
   % creates a res with Super-diffusive particles
   resSupids = res((ismember(res(:,8),supidsLin)),:);
      
   for j = 1:length(supidsLin)
       % index of each element on the matrix
       indfirst(j) = find(resSupids(:,8) == supidsLin(j), 1, 'first') ;
       indlast(j) = find(resSupids(:,8) == supidsLin(j), 1, 'last') ;
       
       lign=[indfirst(:);indlast(:)]; 
       ligns = sort(lign);
       
       tracks = resSupids(ligns,[1:2]);
       tracks = tracks(2:2:end,:)-tracks(1:2:end,:); 
       % tracks length for each particle on x and y
       
       tracksum = (tracks(:,2)).^2+(tracks(:,2).^2); % sum of x and y ^2
       tracksqrt = sqrt(tracksum); % square roots
       tracksqrt = tracksqrt*106 ; % pixel size
   end
   
   % measure of time (nb frames / step)
   [resRec,ia,ic] = unique(resSupids (:,8), 'rows','stable') ; 
   % find a unique value of the particle
   nb_frames = accumarray(ic,1) ; 
   % accumulate the number of frames
   timeSup = nb_frames * MovieTimeStep ; 
   
   % measure the mean speed of formins by embryo
   forminsp = (rdivide(tracksqrt,timeSup))./(1/MovieTimeStep) ; % nm/s 
   mean_forminsp = (mean(forminsp)) ;
 
   save(['mean_forminsp', num2str(i), '.mat'], 'mean_forminsp') ;
   resform = [timeSup, tracksqrt ,forminsp];
   save(['resform', num2str(i) ,'.mat'],'resform');
   
   clearvars -except alphaSub alphaSup barrCc barrI barrRg featsize field goodenough IdivRg maxdisp memory MovieTimeStep num_of_embryos
   
   rmdir KilfoilWorkDir s ; 
 end
 
clear


