%classdef reco_parameter < handle
classdef reco_parameter < matlab.mixin.Copyable
    
    %% Properties
    properties
        scan_name;
        scan_path;
        energy = [];
        sample_detector_distance = [];
        eff_pixel_size = [];
        eff_pixel_size_binned;
        pixel_scaling = 1;
%        shape = [10 10 10];
%         rotation_axis_position
%         rotation_axis_offset
        read_flatcor = 0;
        read_flatcor_path = '';
        read_flatcor_trafo = @(im) im;
        read_sino = 0;
        read_sino_folder = '';
        read_sino_trafo = @(x) (x);
        %%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        raw_roi = [];
        raw_bin = 4;
        im_trafo = '' ;
        proj_range = 1;
        ref_range = 1;
        pixel_filter_threshold_dark = [0.01 0.005];
        pixel_filter_threshold_flat = [0.01 0.005];
        pixel_filter_threshold_proj = [0.01 0.005];
        pixel_filter_radius = [3 3];
        ring_current_normalization = 0;
        image_correlation_method = 'entropy';
        image_correlation_force_calc = 0;
        image_correlation_num_flats = 3;
        image_correlation_area_width = [1 100];
        image_correlation_area_height = [0.2 0.8];
        ring_filter = 0;
        ring_filter_apply_before_stitching = 0;
        ring_filter_method = 'wavelet-fft';
        ring_filter_waveletfft_dec_levels = 2:5;
        ring_filter_waveletfft_wname = 'db25';
        ring_filter_waveletfft_sigma = 2.4;
        ring_filter_jm_median_width = 11;
        strong_absorption_threshold = 1;
        crop_at_rot_axis = 0;
        stitch_projections = 0;
        stitch_method = 'sine';
        %%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phase_retrieval = 0;
        phase_retrieval_apply_before = 1;
        phase_retrieval_post_binning_factor = 1;
        phase_retrieval_method = 'dpc';
        phase_retrieval_reg_par = 1.7;
        phase_retrieval_bin_filt = 0.1;
        phase_retrieval_cutoff_frequ = 2 * pi;
        phase_retrieval_padding = 1;
        dpc_steps = 5;
        dpc_bin = 8;
        %%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tomo = 1;
        tomo_interactive_mode = 1;
        reco_mode =  '3D';
        vol_size = [];
        vol_shape = []
        rot_angle_full_range = [];
        rot_angle_offset = pi;
        angles;
        angle_scaling = 1;
        rot_axis_offset =  0;
        rot_axis_position = [];
        rot_axis_offset_shift;
        tilt_camera = 0;
        tilt_lamino = 0;
        rot_axis_corr_area1 = [];
        rot_axis_corr_area2 = [];
        fbp_filter_type = 'Ram-Lak';
        fbp_filter_freq_cutoff = 1;
        fbp_filter_padding = 1;
        fbp_filter_padding_method = 'symmetric';
        butterworth_filter = 0;
        butterworth_filter_order = 1;
        butterworth_filter_frequ_cutoff = 0.9;
        astra_pixel_size = 1;
        take_neg_log = [];
        algorithm = 'fbp';
        iterations = 40;
        sirt_MinConstraint = [];
        sirt_MaxConstraint = [];
        %%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        write_path = '';
        write_parpath;
        write_to_scratch = 1;
        write_deleteFiles = 0;
        write_beamtimeID_regexp = '';
        write_scan_name_appendix = '';
        write_parfolder = '';
        write_subfolder_flatcor = '';
        write_subfolder_phase_map = '';
        write_subfolder_sino = '';
        write_subfolder_reco = '';
        write_flatcor = 0;
        write_phase_map = 0;
        write_phase_map_path;
        write_phase_appendix;
        write_sino = 0;
        write_sino_path;
        write_sino_phase = 0;
        write_sino_phase_path;
        write_reco = 1;
        write_float = 1;
        write_uint16 = 0;
        write_uint8 = 0;
        write_float_binned = 0;
        write_uint16_binned = 0;
        write_uint8_binned = 0;
        write_reco_binning_factor = 2;
        write_compression_method = 'outlier';
        write_compression_parameter = [0.02 0.02];
        write_uint8_segmented = 0;
        write_reco_path;
        write_figures_path;
        %%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        visual_output = 1;
        interactive_mode_rot_axis_pos = 1;
        interactive_mode_rot_axis_pos_default_search_range = -4:0.5:4; % in binned pixels
        interactive_mode_rot_axis_tilt = 0;
        interactive_mode_rot_axis_tilt_default_search_range = -0.005:0.001:0.005 % radian
        interactive_mode_lamino = 0;
        interactive_mode_angles = 0;
        interactive_mode_angle_scaling_default_search_range = [];
        interactive_mode_slice_number = 0.5;
        interactive_mode_phase_retrieval = 1;
        interactive_mode_phase_retrieval_default_search_range = [];
        %%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        use_cluster = 0;
        poolsize = 0.75;
        poolsize_gpu_limit_factor = 0.7;
        astra_link_data = 1;
        astra_gpu_index = [];
        
        % DYNAMIC
        quick_switch = 0;    
        scan_position = 0;
        offset_shift = [];
        offset_shift_x0;
        offset_shift_x1;        
        vert_shift;
        number_of_stds;
        offset;
        tilt;
        fixed_tilt;
        slice;
        write_reco_phase_path;
        verbose = 0;
        
    end
    
    %% Methods
    methods
        
%         % Rotation axis Offset
%         function obj = set_rotation_axis_offset( obj, value )
%             obj.rotation_axis_offset = value;
%             obj.rotation_axis_position = obj.shape(2) / 2 + value; 
%         end
%         
%         % Rotation axis Position
%        function obj = set_rotation_axis_position( obj, value )
%             obj.rotation_axis_position = value;
%             obj.rotation_axis_offset = value - obj.shape(2) / 2;
%        end
%        
  
%        % Update current object from reference object
%        function obj = update_with( obj, ref )
%            
%            prop_cell = properties( ref );
%            
%            % Loop over properties
%            for nn = 1:numel( prop_cell )
%                prop_name = prop_cell{nn};
%                if ~isempty( ref.(prop_name ) )
%                    obj.(prop_name) = ref.(prop_name);
%                end
%            end
%        end
       
       % Set empty properties with defaults
       function obj = set_empty_prop_with_defaul( obj, ref )
           
           prop_cell = properties( ref );
           
           % Loop over properties
           for nn = 1:numel( prop_cell )
               prop_name = prop_cell{nn};
               if isempty( obj.(prop_name ) )
                   obj.(prop_name) = ref.(prop_name);
               end
           end
       end
       
       
    end
  
    
end
