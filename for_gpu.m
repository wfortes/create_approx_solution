%%%%  

    %% Create the sinogram
    proj_geom = astra_create_proj_geom('parallel',1.0,detector_size,linspace2(0,pi,angle_count));
    vol_geom = astra_create_vol_geom(256,256);

    sinogram_id = astra_mex_data2d( 'store', sinogram ); % Or something like this...

    %% Create the reconstruction
    recon_id = astra_mex_data2d('create','-vol',vol_geom, initial_solution);

    cfg = astra_struct('CGLS_CUDA');
    cfg.ProjectionGeometry = proj_geom;
    cfg.ReconstructionGeometry = vol_geom;
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;

    alg_id = astra_mex_algorithm('create', cfg);
	astra_mex_algorithm('iterate', alg_id, 100);
    reconstruction = astra_mex_data2d('get', recon_id);
    
    astra_mex_data2d('delete',sinogram_id);
    astra_mex_data2d('delete',recon_id);
    astra_mex_algorithm('delete',alg_id);
    disp('created reconstruction');
%%%%%