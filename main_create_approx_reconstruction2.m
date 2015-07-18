
clear all
img_sz_set = [128];%,64,128];%,256,512];
img_index_set = 7;%[1,5,7];
type_set = 1;%[0,1];
[dir_a,dir_b]=mkdirvecs(20);
tol_set = 0;%[0,1,2];
%--------------------------------------------------------------------------
for img_sz = img_sz_set
    if img_sz==8
        N_proj_set = [2,3,4,5,6];
    elseif img_sz==32
        N_proj_set = [2,4,6,8,10,12,14,16];
    elseif img_sz==64
        N_proj_set = [2,4,8,12,16,20,24,28,32];
    elseif img_sz==128
        N_proj_set = 24;%[4,16,24];%[4,8,16,20,24];%,28,32,40,48,56,64];
    elseif img_sz==256
        N_proj_set = [8,16,32,40,48,56,64,72,80,88,96,104];
    elseif img_sz==512
        N_proj_set = [8,16,32,48,64,72,80,88,96,104,112,120,136,152,168,184,200];
    end
    
    for img_index = img_index_set;
        if img_index == 7
            d_set = [0;0.50196;1];
        else
            d_set = [0;1];
        end
        for type_code = type_set;
            for N_proj = N_proj_set;
                for tole = tol_set;
                    if tole == 0
                        tol = 'tol-1';
                    elseif tole == 1
                        tol = 'tol-1:2'; 
                    elseif tole == 2
                        tol = 'tol-09'; 
                    end
%--------------------------------------------------------------------------                
                if type_code == 0
                    type = 'grid';
                    M = mkmatrix(img_sz,img_sz,dir_a(1:N_proj),dir_b(1:N_proj));
                elseif type_code == 1
                    type = 'strip';
                    address = '/export/scratch1/fortes/PhD_files/Load/angles_eq_distr/';
                    M = loadmatrix(address,img_sz,N_proj,type,'matrix');
                end
%
                img = num2str(img_index);
                sz = num2str(img_sz);
                proj = num2str(N_proj);
                
                address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/solind_set/ART/';
                filename = strcat(address,tol,'/',sz,'/solind_set','Im',img,'-sz',sz,'-proj',proj,'-',type);
                load(filename);
                
                thresh = 1/sqrt(img_sz);
%--------------------------------------------------------------------------
                reconstruction = create_approx_solution_piecewise3(M,x,N_proj,d_set,'divide',thresh);
%--------------------------------------------------------------------------
                filename = strcat('d_sol-Im',img,'-sz',sz,'-proj',proj,'-',type,tol,'-divided-t');
                save(filename,'reconstruction');
                end
            end
        end
    end
end