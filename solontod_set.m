clear all
img_sz_set = [64,128,256,512];
img_index_set = 7;%[1,2,3];
type_set = 1;%[0,1];
[dir_a,dir_b]=mkdirvecs(20);
% d_set = [0;1];

for img_sz = img_sz_set
    if img_sz==8
        N_proj_set = [2,3,4,5,6];
    elseif img_sz==32
        N_proj_set = [2,4,6,8,10,12,14,16];
    elseif img_sz==64
        N_proj_set = [2,4,8,12,16,20,24,28,32];
    elseif img_sz==128
        N_proj_set = [4,8,16,20,24,28,32,40,48,56,64];
    elseif img_sz==256
        N_proj_set = [8,16,32,40,48,56,64,72,80,88,96,104];
    elseif img_sz==512
        N_proj_set = [8,16,32,48,64,72,80,88,96,104,112,120,136,152,168,184,200];
    end
    
    for img_index = img_index_set;
        for type_code = type_set;
            for N_proj = N_proj_set;
                
                if type_code == 0
                    type = 'grid';
                    M = mkmatrix(img_sz,img_sz,dir_a(1:N_proj),dir_b(1:N_proj));
                elseif type_code == 1
                    type = 'strip';
                    address = '/export/scratch1/fortes/PhD_files/Load/angles_eq_distr/';
                    M = loadmatrix(address,img_sz,N_proj,type,'matrix');
                end
                
                P = img_read(img_index,img_sz);
                P = reshape(P,img_sz^2,1);
                P = double(P);
                P = P/norm(P,inf); % only for binary images
                Q = M*P;
                
                x = box_constraint_solver(M,Q,1e+5*length(Q),0.9,'ART');
                
                img = num2str(img_index);
                sz = num2str(img_sz);
                proj = num2str(N_proj);
                
                address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/solind_set/ART/tol-09/';
                filename = strcat(address,'solind_set','Im',img,'-sz',sz,'-proj',proj,'-',type);
                save(filename,'x');
                
                N_proj
            end
        end
    end
end

% %%
% 
% clear all
% img_sz_set = [32,64,128,256,512];
% img_index_set = [1,2,3,5];
% type_set = [0,1];
% [dir_a,dir_b]=mkdirvecs(20);
% % d_set = [0;1];
% 
% for img_sz = img_sz_set
%     if img_sz==32
%         N_proj_set = [2,4,6,8,10,12,14,16];
%     elseif img_sz==64
%         N_proj_set = [2,4,8,12,16,20,24,28,32];
%     elseif img_sz==128
%         N_proj_set = [4,8,16,20,24,28,32,40,48,56,64];
%     elseif img_sz==256
%         N_proj_set = [8,16,32,40,48,56,64,72,80,88,96,104];
%     elseif img_sz==512
%         N_proj_set = [8,16,32,48,64,72,80,88,96,104,112,120,136,152,168,184,200];
%     end
%     
%     for img_index = img_index_set;
%         for type_code = type_set;
%             for N_proj = N_proj_set;
%                 
%                 if type_code == 0
%                     type = 'grid';
%                     M = mkmatrix(img_sz,img_sz,dir_a(1:N_proj),dir_b(1:N_proj));
%                 elseif type_code == 1
%                     type = 'strip';
%                     address = '/export/scratch1/fortes/PhD_files/Load/angles_eq_distr/';
%                     M = loadmatrix(address,img_sz,N_proj,type,'matrix');
%                 end
%                 
%                 P = img_read(img_index,img_sz);
%                 P = reshape(P,img_sz^2,1);
%                 P = double(P);
%                 P = P/norm(P,inf); % only for binary images
%                 Q = M*P;
%                 
%                 img = num2str(img_index);
%                 sz = num2str(img_sz);
%                 proj = num2str(N_proj);
%                 
%                 address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/solind_set/ART/tol-1/';
%                 filename = strcat(address,sz,'/solind_set','Im',img,'-sz',sz,'-proj',proj,'-',type);
%                 load(filename);
%                 
%                 %                 N_proj
%                 %                 norm(M*x-Q,inf)
%                 % if isempty(x)
%                 %     fprintf('Im%d-sz%d-proj%d-%s\n',img_index,img_sz,N_proj,type)
%                 % %     error('x nulo')
%                 % end
%                 if min(x)<0 || max(x)>1 || norm(M*x-Q,inf)>0.1 || sum(isnan(x))>0
% %                     for i = 1:length(x)
% %                         if x(i) > d_set(end)
% %                             x(i) = d_set(end);
% %                         elseif x(i) < d_set(1)
% %                             x(i) = d_set(1);
% %                         end
% %                     end
%                     norm(M*x-Q,inf)
%                     min(x)
%                     max(x)
%                     fprintf('Im%d-sz%d-proj%d-%s\n',img_index,img_sz,N_proj,type)
%                 end
%             end
%         end
%     end
% end