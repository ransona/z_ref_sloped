function zref_mid_frame(animalID,ref_name,ch)
if ~exist('animalID')
    error('No animal ID')
end
reg_only = false;
fileID = '';
ch_active = 2; % change to hSI.hChannels.channelsActive once working
%ch = 1; % 1 = green, 2 = red
% prompt for desired spacing of planes in reference stack
plane_spacing = 2; %um
total_planes = 41; %41;
total_range = plane_spacing * (total_planes-1);
% prompt for desired frames per depth
frames_per_depth = 30;

global hSI;
% global hSICtl;
if ~exist('reg_only')
    reg_only = false;
    % 30.05.23 AR reg_only is to be used if you are only registering an
    % already acquired volume
    % fileID = 'Apr-2023 16_53_27';
end

try
    if reg_only == false
        % filename root
        stack_filename_stem = datestr(now);
        stack_filename_stem = strrep(stack_filename_stem,':','_');
        % where to save final tiff
        expDir = fullfile('V:\Local_Repository',animalID,'refz',[ref_name,stack_filename_stem,'.tif']);
        % where to save temporary tifs
        
        % stack_filename_stem = '13-Apr-2023 14_34_52';
        
        expDirTemp = fullfile('V:\Local_Repository',animalID,'refz',stack_filename_stem);
        [ ~, ~ ] = mkdir(fullfile('V:\Local_Repository',animalID,'refz'));
        [ ~, ~ ] = mkdir(expDirTemp);
        x = tic;
        
        logInitialState = hSI.hChannels.loggingEnable;
        hSI.hMotionManager.enable=false;
        hSI.hChannels.loggingEnable=true;
        % start with the fast z profile open that you want to use with the scope
        % in the zoom / zposition that you want to use
        % channel to use
        
        
        hSI.hScan2D.logFilePath = expDirTemp;
        
        % position fast z actuator at zero
        % hSI.hFastZ.positionAbsolute = 0;
        % zeros all motor positions
        hSI.hMotors.setRelativeZero; % this line may change with SI version
        % get z spacing of fast z planes
        fast_z_step = hSI.hStackManager.stackZStepSize;
        fast_z_slices = hSI.hStackManager.numSlices;
        % calculate reference stack z positions
        num_ref_planes = length(0:plane_spacing:fast_z_step);
        ref_stack_z = -(total_range/2):plane_spacing:(total_range/2);
        % set SI volumes to frames_per_depth
        numVolumes_original = hSI.hStackManager.numVolumes;
        hSI.hStackManager.numVolumes = frames_per_depth;
        all_filenames = [];
        for iDepth = 1:length(ref_stack_z)
            disp(['Acquiring volume ',num2str(iDepth),' of ',num2str(length(ref_stack_z))]);
            % set z position
            hSI.hMotors.moveSample([0 0 ref_stack_z(iDepth)]);
            
            % set filename
            % hSI.hScan2D.logFileStem = ['z_ref_',stack_filename_stem,'_',sprintf('%06d', ref_stack_z(iDepth)),'_'];
            hSI.hScan2D.logFileStem = ['z_ref_',stack_filename_stem,'_',sprintf('%06d', iDepth),'_'];
            
            hSI.startGrab;
            % drawnow;
            % pause(3);
            % wait for grab to be complete
            while strcmp(hSI.acqState,'grab')
                drawnow();
            end
        end
        
    else
        % can probably remove this condition:
        % filename root
        stack_filename_stem = fileID;
        % where to save final tiff
        expDir = fullfile('V:\Local_Repository',animalID,'refz',[ref_name,stack_filename_stem,'.tif']);
        expDirTemp = fullfile('V:\Local_Repository',animalID,'refz',stack_filename_stem);
        fast_z_slices = 5;
        fast_z_step = 50;
        
    end
    
    % register
    % load the complete tif of one 5 level stack
    imageFullFileName = dir(fullfile(expDirTemp,['*',stack_filename_stem,'*']));
    [~,idx] = sort([imageFullFileName.datenum]);
    imageFullFileName = imageFullFileName(idx);
    imageFullFileName = {imageFullFileName.name};
    % allocate space to store the registered frames from each depth and the
    % depth value
    all_registered_depths_frames = [];
    all_registered_depths = [];
    total_depths = 0;
    for iFile = 1:length(imageFullFileName)
        disp(['Processing file ',num2str(iFile),'/',num2str(length(imageFullFileName))]);
        info = imfinfo(fullfile(expDirTemp,imageFullFileName{iFile}));
        numberOfPages = length(info);
        all_pages = [];
        current_page = 0;
        % load all pages of tif
        for page_num = ch:ch_active:numberOfPages
            current_page = current_page + 1;
            % Read the kth image in this multipage tiff file.
            all_pages(:,:,current_page) = imread(fullfile(expDirTemp,imageFullFileName{iFile}), page_num);
            % Now process thisPage somehow...
        end
        
        % register pages from same depth of fast z volume
        for iDepth = 3 % only do the middle depth
            disp(['Processing depth ',num2str(iDepth),'/',num2str(fast_z_slices)]);
            slice_frames = all_pages(:,:,iDepth:fast_z_slices:end);
            % test reg function
            % frame1 = double(imread('coins.png'));
            % allfr = cat(3,frame1,circshift(frame1,20,1));
            % regMovie = mean(rapidRegNonPar(allfr,squeeze(allfr(:,:,1))),3);
            % make ref for registration by registering 5 frames - start
            % with 2nd to avoid bad first frame
            slice_frames_ref = mean(rapidRegNonPar(slice_frames(:,:,round(linspace(2,size(slice_frames,3),5))),slice_frames(:,:,round(size(slice_frames,3)/2))),3);
            % register the rest of the frames
            slice_frames_reg = rapidRegNonPar(slice_frames,slice_frames_ref);
            % average the registered frames and allocate them to their depth
            total_depths = total_depths +1;
            all_registered_depths_frames(:,:,total_depths) = mean(slice_frames_reg(:,:,2:end),3);
            % calc depth using fast z slice spacing and file number
            all_registered_depths(total_depths) = ((iDepth - 1)*fast_z_step)+((iFile-1)*plane_spacing);
        end
    end
    
    %     figure;
    %     for iDepth = [1 6 11 16 21 2 7 3 8 4 9 5 10]
    %         imagesc(all_registered_depths_frames(:,:,iDepth))
    %         drawnow
    %         pause(0.3)
    %     end
    
    % sort registered depth frame order depth
    [all_registered_depths_sorted,i] = sort(all_registered_depths);
    all_registered_depths_frames_sorted = all_registered_depths_frames(:,:,i);
    
    % register each depth to the one above it and z score
    % to avoid stange alignment because of the top sections being noise
    % align inside out:
    % from middle to bottom
    for iDepth = floor(size(all_registered_depths_frames_sorted,3)/2):size(all_registered_depths_frames_sorted,3)
        if iDepth > floor(size(all_registered_depths_frames_sorted,3)/2)
            all_registered_depths_frames_sorted(:,:,iDepth) = rapidRegNonPar(all_registered_depths_frames_sorted(:,:,iDepth),all_registered_depths_frames_sorted(:,:,iDepth-1));
        end
        depth_frame = all_registered_depths_frames_sorted(:,:,iDepth);
        depth_frame = depth_frame - min(depth_frame(:));
        depth_frame = depth_frame / max(depth_frame(:));
        depth_frame = depth_frame * 255;
        all_registered_depths_frames_sorted(:,:,iDepth) = depth_frame;
    end
    % from middle to the top
    for iDepth = 1:floor(size(all_registered_depths_frames_sorted,3)/2)
        if iDepth > 1
            all_registered_depths_frames_sorted(:,:,iDepth) = rapidRegNonPar(all_registered_depths_frames_sorted(:,:,iDepth),all_registered_depths_frames_sorted(:,:,iDepth-1));
        end
        depth_frame = all_registered_depths_frames_sorted(:,:,iDepth);
        depth_frame = depth_frame - min(depth_frame(:));
        depth_frame = depth_frame / max(depth_frame(:));
        depth_frame = depth_frame * 255;
        all_registered_depths_frames_sorted(:,:,iDepth) = depth_frame;
    end     
    
%     for iDepth = 1:size(all_registered_depths_frames_sorted,3)
%         if iDepth > 1
%             all_registered_depths_frames_sorted(:,:,iDepth) = rapidRegNonPar(all_registered_depths_frames_sorted(:,:,iDepth),all_registered_depths_frames_sorted(:,:,iDepth-1));
%         end
%         depth_frame = all_registered_depths_frames_sorted(:,:,iDepth);
%         depth_frame = depth_frame - min(depth_frame(:));
%         depth_frame = depth_frame / max(depth_frame(:));
%         depth_frame = depth_frame * 255;
%         all_registered_depths_frames_sorted(:,:,iDepth) = depth_frame;
%     end

    f = figure;
    while true
        for iDepth = 1:size(all_registered_depths_frames_sorted,3)
            imagesc(all_registered_depths_frames_sorted(:,:,iDepth),[0 256])
            colorbar
            drawnow
            pause(0.2)
        end
        if ~strcmp(questdlg('Replay stack?'),'Yes')
            close(f)
            break
        end
    end            
    
    all_registered_depths_frames_sorted = uint16(all_registered_depths_frames_sorted);
    % output as a registration object for SI online motion correction
    imwrite(squeeze(all_registered_depths_frames_sorted(:,:,1)), expDir);
    for iFrame = 2:size(all_registered_depths_frames_sorted,3)
        imwrite(squeeze(all_registered_depths_frames_sorted(:,:,iFrame)), expDir, 'WriteMode', 'append');
    end
    disp(['Time taken = ',num2str(toc(x))]);
catch err1
    disp('Error processing stack')
    disp(err1.identifier)
end

if reg_only == false
    % put settings back to start state
    hSI.hStackManager.numVolumes = 99999999;
    hSI.hChannels.loggingEnable=logInitialState;
    hSI.hScan2D.logFileStem = '';
    hSI.hScan2D.logFilePath = '';
    % set motor to where it was
    hSI.hMotors.moveSample([0 0 0]);
end
