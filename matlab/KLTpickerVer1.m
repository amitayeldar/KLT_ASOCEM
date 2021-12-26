function KLTpickerVer1(micrograph_addr,output_dir,particle_size,num_of_particles,num_of_noise_images,use_ASOCEM,gpu_use)
% 
% KLT picker
% 
%
% Amitay Eldar, Nov 2019

% Input:
% micrograph_addr       address of the micrographs .mrc files.
% output_dir            address to save coordinate file
% particle_size         particle size in pixels
% num_of_particles      number of particles to pick per micrograph. if set to -1
%                       then pick all particles
% num_of_noise_images   number of noise images to pick per micrograph
% gpu_use               if set to 1 then use the GPU.


% Output:
% coordinate files .box in output_dir
% Picking summery text file in output_dir.
         
%% Initielizing parameters
mgScale = 101/particle_size; % scaling to particle size of 101 should be odd number  pix which seems to give good results in practice
files = dir([micrograph_addr,'/','*.mrc']);
numOfMicro = size(files,1);
patchSz = floor(0.8*mgScale*particle_size);% apprx Par size after downSampling
patchSzFun = floor(0.4*mgScale*particle_size); % disk Sz for computing the eigen func
patchSzPickBox = floor(1*mgScale*particle_size); % box size
maxIter = 6*10^4; % max iterations for psd apprx
MaxNumOfFun= 400; % max eigen function to use.
maxOrder= 100;% order of maximum eigen function
precentOfEig = 0.99; % how many eigen function to develope
thresh = 0; % threshold for the picking.
showFig = 0; % show figs if set to 1.
preProcessStage = 1; % 1 to creat variables for all micrographs
microNames = cell(numOfMicro,1);
pickedParPerMic = zeros(numOfMicro,1);
pickedNoisePerMic = zeros(numOfMicro,1);

% automatic part
if mod(patchSz,2) == 0
    patchSz = patchSz -1;
end
if mod(patchSzFun,2) == 0
    patchSzFun = patchSzFun -1;
end

coordinatsPathParticle =[output_dir,'pickedParticles','ParticleSize',num2str(particle_size)];
coordinatsPathNoise =[output_dir,'pickedNoise','ParticleSize',num2str(particle_size)];
if ~exist(coordinatsPathParticle, 'dir')
   mkdir(coordinatsPathParticle)
end
if ~exist(coordinatsPathNoise, 'dir')
    if num_of_noise_images~=0
        mkdir(coordinatsPathNoise)
    end
end


%% PreProcess
disp("Starting preprocess");
if preProcessStage == 1
    T = 1; % Nyquist sampling rate
    bandLimit = pi/T;
    radMax =  floor((patchSzFun-1)/2);
    x = double(-radMax:1:radMax);
    [X,Y] = meshgrid(x);
    radMat = sqrt(X.^2+Y.^2); theta = atan2(Y,X);
    szRadMat = size(radMat);
    name=[' max order ',num2str(maxOrder),' precentOfEig ',num2str(precentOfEig)]; % for figures
    [rSamp,~,idxRsamp] = uniquetol(radMat(:),1e-14);
    theta = theta(:);
    numOfQuadKer =2^7;
    numOfQuadNys = 2^7;
    [rho,quadKer] = lgwt(numOfQuadKer,0,bandLimit);
    rho = flipud(double(rho)); quadKer = flipud(double(quadKer)); 
    [r,quadNys] = lgwt(numOfQuadNys,0,radMax);
    r = flipud(double(r)); quadNys = flipud(double(quadNys));
    % mat for finding eig
    rr = r*r';
    rRho = r*rho'; 
    % mat for sampeling eigfunction
    One = double(ones(length(rSamp),1));
    rSampr = One*r';
    rSampRho = rSamp*rho';
    sqrtrSampr = sqrt(rSampr);
    JrRho = double(zeros(size(rRho,1),size(rRho,2),length(0:maxOrder-1)));
    Jsamp = double(zeros(size(rSampRho,1),size(rSampRho,2),length(0:maxOrder-1)));
    cosine = double(zeros(length(theta),length(0:maxOrder-1)));
    sine = double(zeros(length(theta),length(0:maxOrder-1))); 
    parfor N=0:maxOrder-1
       JrRho(:,:,N+1) = double(besselj(N,rRho)); 
       Jsamp(:,:,N+1) = double(besselj(N,rSampRho));
        if N~=0
            cosine(:,N+1) = double(cos(N*theta));
            sine(:,N+1) = double(sin(N*theta));
        end
    end
    
end
disp('Preprocess finished');

%% main section
disp("Starting particle picking from micrographs")

% Initialize queue for progress messages
print_progress % Call without paramters to reset all progress variables
progressQ=parallel.pool.DataQueue;
afterEach(progressQ, @print_progress); % function defined at the end

parfor expNum = 1:numOfMicro
    startT=clock;
    [~, microName] = fileparts(files(expNum).name);
    mgBig = ReadMRC([files(expNum).folder,'/',files(expNum).name]);
    mgBig = double(mgBig);
    mgBig = rot90(mgBig);
    mgBigSz = size(mgBig);
    mg = cryo_downsample(mgBig,[floor(mgScale*size(mgBig,1)),floor(mgScale*size(mgBig,2))]);
    if mod(size(mg,1),2) == 0 % we need odd size
        mg = mg(1:end-1,1:end);
    end
    if mod(size(mg,2),2) == 0 % we need odd size
        mg = mg(1:end,1:end-1);
    end
 
    
    %% contamination removal using ASOCEM

    if use_ASOCEM==1
        I0 = imgaussfilt(mgBig,1);
        maxIterAsocem=300;
        downscale_size = 600;
        contamination_criterion = 0; % that means by size.
        fast_flag = 2; % that means fast.
        area_size = 5; % after down scaling to 400*400
        % run ASOCEM
        stop_d =0;
        while stop_d == 0
            if rem(downscale_size,area_size) == 0
                if rem(downscale_size,2) ==1
                   stop_d = 1;
                end
            end
            if stop_d == 0
                downscale_size = downscale_size -1;
            end
        end

        [phi] = ASOCEM(I0,downscale_size,area_size,contamination_criterion,fast_flag,maxIterAsocem);
        if phi==ones(size(phi)) % all is contamination dont pick
             continue
        end
        % get rid of the edges
        scalingSz = downscale_size/max(size(I0));
        d = max(3,ceil(scalingSz*particle_size/8));
        phi(1:d,:) = 0; phi(end-d+1:end,:)=0; phi(:,1:d)=0; phi(:,end-d+1:end)=0;
        % get rid of blubs the size of the particle
        phi_seg = imbinarize(zeros(size(phi)));%     min_bulb_size = floor(2*scalingSz*particle_size/2);
        se_erod = strel('square',max(area_size,ceil(scalingSz*particle_size/6)));
        phi_erod = imerode(phi>0,se_erod);
        CC = bwconncomp(phi_erod,8);
        for i =1:size(CC.PixelIdxList,2)
            if size(CC.PixelIdxList{i},1)> (scalingSz*particle_size)^2
                    phi_seg(CC.PixelIdxList{i})=1;
            end
        end 
        phi_seg = imdilate(phi_seg,se_erod);
        f=figure('visible', 'off');
        subplot(1,2,1); imshow(cryo_downsample(I0,200),[]);
        subplot(1,2,2); imshow(imresize(phi_seg,[200,200]),[]);
        mkdir([output_dir,'/AsocamFigs']);
        saveas(f,[output_dir,'/AsocamFigs/',microName,'.jpg'])
%         WriteMRC(imresize(phi_seg,size(mgBig)),1,[output_dir,'/AsocamFigs/',microName,'.mrc'])
        phi_seg = imresize(phi_seg,size(mg));
    else
        phi_seg = zeros(size(mg));
    end
    
    %% normalization
    tmp = mg(phi_seg<=0);
    mg = mg - mean(tmp(:));
    mg = mg/norm(tmp - mean(tmp(:)),'fro'); % normalization; 
    mcSz = size(mg);

    %% Cutoff filter
    noiseMc = mg; 
    bandPass1d = fir1(patchSz-1, [0.05 0.95]);
    bandPass2d  = ftrans2(bandPass1d); %% radial bandpass
    if gpu_use == 1
        noiseMc = imfilter(gpuArray(noiseMc), bandPass2d);
        noiseMc = gather(noiseMc);
    else
        noiseMc = imfilter(noiseMc, bandPass2d);
    end
    %% Estimating particle and noise RPSD
    [apprxCleanPsd,apprxNoisePsd,~,R,~,stopPar] = rpsd_estimation(noiseMc,phi_seg,patchSz,maxIter,gpu_use);
    if stopPar==1 % maxIter happend, skip to next micro graph
        continue
    end
    if showFig==1
        figure('visible','on');
        plot(R*bandLimit,apprxCleanPsd,'LineWidth',4)
        grid on
        figName = ('apprx Clean Psd first stage ');
        legend('apprxCleanPsd')
        title(figName,'FontSize',9);
        figure('visible','on');
        plot(R*bandLimit,apprxNoisePsd,'LineWidth',4)
        grid on
        figName = ('apprx Noise Psd first stage ');
        legend('apprxNoisePsd')
        title(figName,'FontSize',9);
    end
   
    %% PreWhitening the micrograph
    apprxNoisePsd = apprxNoisePsd +(median(apprxNoisePsd)*(10^-1));% we dont want zeros
    [mgPrewhite] = prewhite(mg,apprxNoisePsd,apprxCleanPsd,mcSz,patchSz,R);
    
    
    %% Re estimating particle and noise RPSD
    % normalization
    tmp = mgPrewhite(phi_seg<=0);
    mgPrewhite = mgPrewhite - mean(tmp(:));
    mgPrewhite = mgPrewhite/norm(tmp - mean(tmp(:)),'fro'); % normalization;
    %% Cutoff filter
    noiseMc = mgPrewhite; 
    bandPass1d = fir1(patchSz-1, [0.05 0.95]);
    bandPass2d  = ftrans2(bandPass1d); %% radial bandpass
    if gpu_use == 1
        noiseMc = imfilter(gpuArray(noiseMc), bandPass2d);
        noiseMc = gather(noiseMc);
    else
        noiseMc = imfilter(noiseMc, bandPass2d);
    end
    [apprxCleanPsd,apprxNoisePsd,noiseVar,R,~,stopPar] = rpsd_estimation(noiseMc,phi_seg,patchSz,maxIter,gpu_use);
    if stopPar==1 % maxIter happend, skip to next micro graph
        continue
    end 
    if showFig==1
        figure('visible','on');
        plot(R*bandLimit,apprxCleanPsd,'LineWidth',4)
        grid on
        figName = ('apprx Clean Psd second stage ');
        legend('apprxCleanPsd')
        title(figName,'FontSize',9);
        figure('visible','on');
        plot(R*bandLimit,apprxNoisePsd,'LineWidth',4)
        grid on
        figName = ('apprx Noise Psd second stage ');
        legend('apprxNoisePsd')
        title(figName,'FontSize',9);
    end


    %% Constructing the KLTpicker templates
    psd = abs(triginterp(rho,bandLimit*R,apprxCleanPsd));
    if showFig==1
        figure('visible','on');
        plot(rho,psd);
        figName = ['Clean Sig Samp at nodes ',name];
        title(figName);
    end
    [eigFun,eigVal] = construct_klt_templates(rho,quadKer,quadNys,rr,sqrtrSampr,JrRho,Jsamp,cosine,sine,numOfQuadNys,maxOrder,psd,precentOfEig,idxRsamp,gpu_use);

   if size(eigFun,2) < MaxNumOfFun 
        numOfFun = size(eigFun,2);
   else
        numOfFun = MaxNumOfFun;
   end
%     for i = 1:size(eigFun,2)
%        
%         tmpFun(:,:,i) = reshape(eigFunStat(:,i),patchSzFun,patchSzFun);
%     end


    %% particle detection
    [numOfPickedPar,numOfPickedNoise] = particle_detection(mgPrewhite,phi_seg,eigFun,eigVal,numOfFun,noiseVar,mcSz,mgScale,radMat,mgBigSz,patchSzPickBox,patchSzFun,num_of_particles,num_of_noise_images,coordinatsPathParticle,coordinatsPathNoise,microName,thresh,gpu_use);
    microNames{expNum} = microName;
    pickedParPerMic(expNum) = numOfPickedPar;
    pickedNoisePerMic(expNum) = numOfPickedNoise;
    % Report progress
    endT=clock;
    data=struct;
    data.i=expNum; 
    data.n_mics=numOfMicro; 
    data.t=etime(endT,startT);
    data.numOfPickedPar=numOfPickedPar;
    send(progressQ,data);

end
disp("Finished the picking successfully");
message = ['Picked ', num2str(sum(pickedParPerMic)),' particles and ',num2str(sum(pickedNoisePerMic)),' noise images out of ',num2str(numOfMicro),' micrographs.'];
disp(message);
% creating a text file the summerizes the number of picked particles and noise
disp("Writing Picking Summery at the output path");
pickingSummery = fopen(fullfile(output_dir,['pickingSummery','.txt']),'w');
fprintf(pickingSummery,'%s\n','Picking Summery');
fprintf(pickingSummery,'%s\n',message);
fprintf(pickingSummery,'%s\n','');
fprintf(pickingSummery,'%s\n','Picking per micrograph:');
fprintf(pickingSummery,'%s\n','Micrographs name #1');
fprintf(pickingSummery,'%s\n','Number of picked particles #2');
fprintf(pickingSummery,'%s\n','Number of picked noise images #3');
fprintf(pickingSummery,'%s\n','--------------------------------');
for i = 1:numOfMicro
    fprintf(pickingSummery,'%s\t%i\t%i\n',microNames{i},pickedParPerMic(i),pickedNoisePerMic(i));
end
fclose(pickingSummery);
disp("The KLT picker has finished");
end

function print_progress(data)
    persistent tot_mic
    persistent tot_par
    persistent start_time
    persistent remaining_time
    
    % If no arguments are given to the function, then reset all variables.
    if nargin==0
        tot_mic=[];
        tot_par=[];        
        start_time=[];
        remaining_time=[];
        return
    end        
    
    if isempty(tot_mic)
        tot_mic = 0;  % Total number of micrographs processed
    end
    
    if isempty(tot_par)
        tot_par = 0;  % Total number of particles picked
    end
        
    if isempty(start_time)
        start_time=clock; %Timestamp of starting time
    end

    if isempty(remaining_time)
        remaining_time=0; % Estimated time left for processing
    end
    
    tot_mic=tot_mic+1;
    tot_par=tot_par+data.numOfPickedPar;
    tot_time=etime(clock,start_time);
    avg_time=tot_time/tot_mic; % Average processing time per micrograph
        
    p=gcp('nocreate');
    if tot_mic > 2*p.NumWorkers
        t_est = (data.n_mics-tot_mic)*avg_time; % Current estimate for remaining time.
        remaining_time= 0.6*remaining_time + 0.4 * t_est; % Smooth remaining time.
        fprintf('Done micrograph %04d (%d/%d) in %3.0f secs (ETA %.0f mins). Picked so far %d particles.\n',...
            data.i,tot_mic,data.n_mics,data.t,remaining_time/60,tot_par);
    else % Don't print ETA for the first micrographs
        remaining_time=(data.n_mics-tot_mic)*avg_time; % Not enough 
            % micrographs yet to smooth remaining time, so make a crude
            % estimate.
        fprintf('Done micrograph %04d (%d/%d) in %3.0f secs (ETA [still estimating]). Picked so far %d particles.\n',...
            data.i,tot_mic,data.n_mics,data.t,tot_par);
    end

end
