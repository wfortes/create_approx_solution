
clear all
img_sz_set = [32,128];%,256,512];
img_index_set = [1,5,7];
type_set = 1;%[0,1];
[dir_a,dir_b]=mkdirvecs(20);
tol_set = 0;%[0,1,2];
rounded = 0;%norm(M*(R-x),inf);
%--------------------------------------------------------------------------
for img_sz = img_sz_set
    tau = 1/img_sz;%sqrt(img_sz);
    
    if img_sz==8
        N_proj_set = [2,3,4,5,6];
    elseif img_sz==32
        N_proj_set = [2,4,6,8,10,12,14,16];
    elseif img_sz==64
        N_proj_set = [2,4,8,12,16,20,24,28,32];
    elseif img_sz==128
        N_proj_set = [4,8,16,20,24];%,28,32,40,48,56,64];
    elseif img_sz==256
        N_proj_set = [8,16,32,40,48,56,64,72,80,88,96,104];
    elseif img_sz==512
        N_proj_set = [8,16,32,48,64,72,80,88,96,104,112,120,136,152,168,184,200];
    end
    
    for img_index = img_index_set;
        if img_index == 7
            d_set = [0;0.50196;1];%0.5019607843
            z = 0.50196;
        else
            d_set = [0;1];
            z = 1;
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
                
%                 address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/reconstruction/';
%                 filename = strcat(address,'d_sol-','Im',img,'-sz',sz,'-proj',proj,'-',type);
%                 address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/divided/t/';
                address = '/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/thresh/t_img_sz/';
                filename = strcat(address,'d_sol-','Im',img,'-sz',sz,'-proj',proj,'-',type,tol,'-t2');
                load(filename);
                x=reconstruction;
                
                
                P = img_read(img_index,img_sz);
                P = reshape(P,img_sz^2,1);
                P = double(P);
                P = P/norm(P,inf);
                
                dif = abs(M*(P-x));
%                 [a,b,c] = find(dif > (N_proj-1)*largest_dif);
%--------------------------------------------------------------------------                
%                 if norm(dif,inf) <= sqrt(norm(M,1))*z% +1e-1+ (norm(M,inf)-N_proj)*tau + rounded
%                     fprintf('%d:square root\n',img_index)
%                 elseif norm(dif,inf) <= norm(M,1)*z/2% +1e-1+ (norm(M,inf)-N_proj)*tau + rounded
%                     fprintf('%d:1/2\n',img_index)
%                 elseif norm(dif,inf) <= (norm(M,1)-2)*z% +1e-1+ (norm(M,inf)-N_proj)*tau + rounded
%                     fprintf('%d:-2\n',img_index)
%                 elseif norm(dif,inf) <= (norm(M,1)-1.5)*z%+1e-1 + (norm(M,inf)-N_proj)*tau + rounded
%                     fprintf('%d:-1.5\n',img_index)
%                 elseif norm(dif,inf) <= (norm(M,1)-1)*z% +1e-1+ (norm(M,inf)-N_proj)*tau + rounded
%                     fprintf('%d:-1\n',img_index)
%                 else
                    if norm(dif,inf) <= norm(M,1)*z +1e-1+ (norm(M,inf)-N_proj)*tau + rounded
                    fprintf('%d:just\n',img_index)
                else
                    fprintf('%d:%g\n',img_index,(abs(max(abs(dif)) - (N_proj)*z)+1e-1+ (norm(M,inf)-N_proj)*tau + rounded)/z)
                    %             if norm(data.dif,inf) > N_proj*a*data.z/2+1e-13 %+ (norm(M,inf)-N_proj+1)*tau + rounded
                    %                 fprintf('%d:%g\n',img_index,(abs(max(abs(data.dif)) - (N_proj*a)*data.z)/2)/((N_proj*a)*data.z)/2)
                end
                %             fprintf('%d\n%g\n%g\n%g\n',N_proj,(N_proj-1)*data.z,max(abs(data.dif)),norm(data.R-data.reconstruction,1)/length(data.dif))
%--------------------------------------------------------------------------
                end
            end
        end
    end
end
    % piece=4;
    % n = 64;
    % mm = n/piece^2;
    % m = sqrt(mm);
    % A=[];
    % B=[];
    % x=1:64;
    % for i_pieces_columns = 1:m
    %     for i_pieces_rows = 1:m
    %         Ix = [];
    %         pixels_in_piece = [];
    %
    %         for i_row = 0:piece-1
    %             pixels_in_piece_aux = i_row*piece*m+1+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1):i_row*piece*m+piece+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1);
    %             pixels_in_piece = union(pixels_in_piece,pixels_in_piece_aux);
    %         end
    %
    % %         for i_aux = 1:piece^2
    %
    %                 a = reshape(x(pixels_in_piece),piece,piece)
    %                 A = [A;a];
    % %         end
    %     end
    %      B = [B A];
    %     A=[];
    % end
    %                 B
    %                 reshape(B,piece^2*mm,1)