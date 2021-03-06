function KLTpicker_start
%
% KLTpicker_start   Start KLT picker processing workflow
%
% KLTpicker_start
% Gather all information required to start the particle picking process out the micrographs.
% This is the first command to be called in any processing workflw.
% 
% Amitay Eldar, November 2019.

if ~isdeployed % Only run in a MATLAB session
    [basedir,~,~]=fileparts(mfilename('fullpath'));
    addpath(fullfile(basedir,'matlab')); % set up MATLAB path
end

micrograph_addr='';
while isempty(micrograph_addr)
    micrograph_addr =fmtinput('Enter full path of micrographs MRC file: ','','%s');
%     if exist(micrograph_addr,'file')~=7
    if isempty(dir([micrograph_addr,'/*.mrc']))
        fprintf('MRC file does not exist.\n');
        micrograph_addr='';
    end
end

output_dir =fmtinput('Enter full path of output directory: ','','%s');
if ~strcmp(output_dir(end),'/')
    output_dir = [output_dir,'/'];
end
    
if ~exist(output_dir,'dir') % Do we need to create directory?
    message='Output directory does not exist. Create?';
    do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    if do_create==1
        mkdir(output_dir);
    end
end

particle_size='';
while isempty(particle_size)
    particle_size_str =fmtinput('Enter the particle size in pixels: ','','%s');
    particle_size = str2double(particle_size_str);
    if mod(particle_size,1)~=0
        fprintf('particle size should be a natural number.\n');
        particle_size='';
    end
    if particle_size<0
        fprintf('particle size should be a natural number.\n');
        particle_size='';
    end
end


num_of_particles='';
message='pick all particles?';
do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_create==1
    num_of_particles = -1;
end
while isempty(num_of_particles)
   num_of_particles_str =fmtinput('How many particles to pick: ','','%s');
   num_of_particles = str2double(num_of_particles_str);
    if mod(num_of_particles,1)~=0
        fprintf('number of particles to pick should be a natural number.\n');
        num_of_particles='';
    end
    if num_of_particles<0
        fprintf('number of particles to pick should be a natural number.\n');
        num_of_particles='';
    end
end

num_of_noise_images='';
message='Pick noise images?';
do_create=multichoice_question(message,{'Y','N'},[ 0, 1],'N');
if do_create==1
    num_of_noise_images=0;
end
while isempty(num_of_noise_images)
   num_of_noise_images_str =fmtinput('How many noise images to pick: ','','%s');
   num_of_noise_images = str2double(num_of_noise_images_str);
    if mod(num_of_noise_images,1)~=0
        fprintf('number of noise images to pick should be a natural number.\n');
        num_of_noise_images='';
    end
    if num_of_noise_images<0
        fprintf('number of noise images to pick should be a natural number.\n');
        num_of_noise_images='';
    end
end

message='Do you want to use ASOCEM for contamination removal?';
do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_create==1
    addpath('./matlab/ASOCEM_functions');
    use_ASOCEM = 1;
    smoothing_term='';
while isempty(smoothing_term)
    smoothing_term_str =fmtinput('Enter the smoothing_term 1 for max smoothing and 0 not to smooth: ','','%s');
    smoothing_term = str2double(smoothing_term_str);
    if or(smoothing_term<0,smoothing_term>1)
        fprintf('smoothing_term should be a number between 0 to 1.\n');
        smoothing_term='';
    end
end
else
    use_ASOCEM = 0;
    smoothing_term = 0;
end

message='Do you want to use the GPU?';
do_create=multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
if do_create==1
    gpu_use = 1;
else
    gpu_use = 0;
end


KLTpickerVer1(micrograph_addr,output_dir,particle_size,num_of_particles,num_of_noise_images,use_ASOCEM,smoothing_term,gpu_use)

