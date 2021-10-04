
function ASOCEM_ver1(micrograph_addr,output_dir,particle_size,downscale_size,area_size,contamination_criterion,fast_flag)
% parameters
files = dir([micrograph_addr,'/','*.mrc']);
numOfMicro = size(files,1);
maxIterAsocem=600;
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
for expNum = 1:numOfMicro
    [~, microName] = fileparts(files(expNum).name);
    mgBig = ReadMRC([files(expNum).folder,'/',files(expNum).name]);
    mgBig = double(mgBig);
%     mgBig = rot90(mgBig);
    I0 = mgBig;
    % run ASOCEM
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
        if size(CC.PixelIdxList{i},1)> (scalingSz*2*particle_size)^2 
                phi_seg(CC.PixelIdxList{i})=1;
        end
    end 
    phi_seg = imdilate(phi_seg,se_erod);
    % resizing phi_seg to original image
    %% we want to save phi_seg as a binary mat file
    f=figure('visible', 'off');
    subplot(1,2,1); imshow(cryo_downsample(I0,size(phi_seg)),[]);
    subplot(1,2,2); imshow(phi_seg,[]);
    saveas(f,[output_dir,microName,' _imageSz_',num2str(downscale_size),' area_ ',num2str(area_size),' fast_flag ',num2str(fast_flag),'.jpg'])
    phi_seg = imresize(phi_seg,size(mgBig));

end
end