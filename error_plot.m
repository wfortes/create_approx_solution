
clear all
img_sz_set = [32,128];%,64,128];%,256,512];
img_index_set = [1,5,7];
type_set = 1;%[0,1];
[dir_a,dir_b]=mkdirvecs(20);
tol_set = 0;%[0,1,2];
method_set = 2;%[1,2,3,4];
%--------------------------------------------------------------------------
for img_sz = img_sz_set
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
            d_set = [0;0.50196;1];
            d = 0.50196;
        else
            d_set = [0;1];
            d = 1;
        end
        for type_code = type_set;
            for tole = tol_set;
                if tole == 0
                    tol = 'tol-1';
                    ntol = 0.1;
                elseif tole == 1
                    tol = 'tol-1:2';
                    ntol = 0.5;
                elseif tole == 2
                    tol = 'tol-09';
                    ntol = 0.9;
                end
                for method = method_set
                    aux = 0;
                Bound = zeros(length(N_proj_set),1);
                Actual_p_dif = zeros(length(N_proj_set),1);
                Image_error = zeros(length(N_proj_set),1);
                for N_proj = N_proj_set;
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
                    P = img_read(img_index,img_sz);
                    P = reshape(P,img_sz^2,1);
                    P = double(P);
                    P = P/norm(P,inf);
                    %
                    img = num2str(img_index);
                    sz = num2str(img_sz);
                    proj = num2str(N_proj);
                    address2 ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/graphs/';
                    
                    if method == 1
                        address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/reconstruction/';
                        filename = strcat(address,tol,'/d_sol-Im',img,'-sz',sz,'-proj',proj,'-',type,tol);
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',type,tol,'.fig');
                        thresh = 0;
                    elseif method == 2
                        address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/thresh/t3/';
                        filename = strcat(address,'d_sol-Im',img,'-sz',sz,'-proj',proj,'-',type,tol,'-t3');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',type,tol,'-t3.fig');
                        thresh = 1/img_sz;%sqrt(img_sz);
                    elseif method == 3
                        address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/divided/';
                        filename = strcat(address,tol,'/d_sol-Im',img,'-sz',sz,'-proj',proj,'-',type,tol,'-divided');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',type,tol,'-d.fig');
                        thresh = 0;
                    elseif method == 4
                        address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/divided/t/';
                        filename = strcat(address,'d_sol-Im',img,'-sz',sz,'-proj',proj,'-',type,tol,'-divided-t');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',type,tol,'-dt.fig');
                        thresh = 1/sqrt(img_sz);
                    end
                    load(filename);
                    aux =aux + 1;
%--------------------------------------------------------------------------
                    Bound(aux,1) = (d*norm(M,1) + ntol + (norm(M,inf)-norm(M,1))*thresh);%/norm(M*P);
                    Actual_p_dif(aux,1) = norm(M*(P-reconstruction),inf);%/norm(M*P);
                    Image_error(aux,1) = norm(P-reconstruction,2);%/length(P);
%--------------------------------------------------------------------------
                end
                x1 = N_proj_set;
                y0 = Bound;
                y1 = Actual_p_dif;
                y2 = Image_error;
                
                figure
                h1 = line(x1,y0);
                h2 = line(x1,y1);
                ax1 = gca;
                ax2 = axes('Position',get(ax1,'Position'),...
                'XAxisLocation','bottom',...
                'YAxisLocation','right',...
                'Color','none',...
                'XTick',[],...
                'XColor','k','YColor','r');
                h3 = line(x1,y2);
                set(get(ax1,'Ylabel'),'String','Projection distance','fontsize',17.45,'Color','k') 
                set(get(ax2,'Ylabel'),'String','Image distance','fontsize',17.45,'Color','r') 
                set(ax1,'fontsize',15,'ycolor','k')
                set(ax2,'fontsize',15,'ycolor','r')
                set(h1,'Color','k','LineStyle','--','Marker','+','LineWidth',2,'MarkerSize',10);
                set(h2,'Color','k','Marker','x','LineWidth',2,'MarkerSize',10);
                set(h3,'Color','r','Marker','s','LineWidth',2,'MarkerSize',10);

                legend([h1, h2, h3],'B','Pd','Id','location','best')
%                 text(6.5,0.1,'Number of angles','fontsize',17.45)
%                 xlabel('Number of angles','fontsize',17.45)
                
                %
                saveas(h1,filename2);
                clear h1 h2 h3
                close all
                end
            end
        end
    end
end

%%



clear all
img_sz_set = [32,128];%,64,128];%,256,512];
img_index_set = [1,5,7];
type_set = 1;%[0,1];
[dir_a,dir_b]=mkdirvecs(20);
tol_set = 0;%[0,1,2];
method_set = 2;%[1,2,3,4];
%--------------------------------------------------------------------------
for img_sz = img_sz_set
%     if img_sz==8
%         N_proj_set = [2,3,4,5,6];
%     elseif img_sz==32
%         N_proj_set = [2,4,6,8,10,12,14,16];
%     elseif img_sz==64
%         N_proj_set = [2,4,8,12,16,20,24,28,32];
%     elseif img_sz==128
%         N_proj_set = [4,8,16,20,24];%,28,32,40,48,56,64];
%     elseif img_sz==256
%         N_proj_set = [8,16,32,40,48,56,64,72,80,88,96,104];
%     elseif img_sz==512
%         N_proj_set = [8,16,32,48,64,72,80,88,96,104,112,120,136,152,168,184,200];
%     end
    
    for img_index = img_index_set;
        if img_index == 7
            d_set = [0;0.50196;1];
            d = 0.50196;
        else
            d_set = [0;1];
            d = 1;
        end
        for type_code = type_set;
            if type_code == 0
                typ = 'grid';
            elseif type_code == 1
                typ = 'strip';
            end
            for tole = tol_set;
                if tole == 0
                    tol = 'tol-1';
                    ntol = 0.1;
                elseif tole == 1
                    tol = 'tol-1:2';
                    ntol = 0.5;
                elseif tole == 2
                    tol = 'tol-09';
                    ntol = 0.9;
                end
                for method = method_set
%                 for N_proj = N_proj_set;
%--------------------------------------------------------------------------
                    img = num2str(img_index);
                    sz = num2str(img_sz);
                    address ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/graphs/';
                    address2 ='/ufs/fortes/Desktop/PhD_m_files/tomography/create_approx_solution/graphs/pdf/';
                    
                    if method == 1
                        filename = strcat(address,'graph',img,'-sz',sz,'-',typ,tol,'.fig');                        
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',typ,tol);
                        thresh = 0;
                    elseif method == 2
                        filename = strcat(address,'graph',img,'-sz',sz,'-',typ,tol,'-t3.fig');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',typ,tol,'-t3');
                        thresh = 2/sqrt(img_sz);
                    elseif method == 3
                        filename = strcat(address,'graph',img,'-sz',sz,'-',typ,tol,'-d.fig');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',typ,tol,'-d');
                        thresh = 0;
                    elseif method == 4
                        filename = strcat(address,'graph',img,'-sz',sz,'-',typ,tol,'-dt.fig');
                        filename2 = strcat(address2,'graph',img,'-sz',sz,'-',typ,tol,'-dt');
                        thresh = 1/sqrt(img_sz);
                    end
%--------------------------------------------------------------------------
                    open(filename)
                    set(gca,'color','none')
                    set(gcf,'color','none')
                    export_fig(filename2,'-pdf','-painters','-native');
                    close all
%                 end
                end
            end
        end
    end
end