function [TFNBS_Result] = gretna_TFNBS(Data_path, File_filter, M, Mask_net, Pthr, E, H, Cova_path)

%==========================================================================
% This function is used to perform the TFNBS algorithm to search
% connections that show significant group effects for one-way experimental
% design (two or more groups).
%
%
% Syntax: function [TFNBS_Result] = gretna_TFNBS(Data_path, File_filter, M, Mask_net, Pthr, E, H, Cova_path)
%
% Inputs:
%      Data_path:
%                N*1 cell with each cell denoting the directory where
%                connectivity matrices are sorted for each group.
%    File_filter:
%                N*1 cell with each cell denoting the prefix of those
%                connectivity matrices to be processed for each group.
%              M:
%                The number of permutations.
%       Mask_net:
%                The matrix mask containing 0 and 1 and only connections
%                corresponding to 1 are fed in the TFNBS computation.
%           Pthr:
%                The FWE-corrected p-value threshold.
%              E:
%                Extension enhancement parameters(default = 0.5).
%              H:
%                Height enhancement parameters(default = 2).
% Cova_path (opt):
%                N*1 cell with each cell denoting the directory and filename
%                of covariates for each group. NB, the group order should
%                be the same to that of Data_path.
%
% Output:
%   TFNBS_Result.:
%           Fmat:
%                F values for each edge of interest.
%         P_Fmat:
%                P values of Fmat.
%        DF_Fmat:
%                Degree of freedom of Fmat.
%           Tmat:
%                Post hoc T values for each edge of interest. For n groups,
%                the order is: g(1)-g(2), g(1)-g(3), ..., g(1)-g(n),
%                g(2)-g(3), ..., g(2)-g(n), ..., g(n-1)-g(n).
%         P_Tmat:
%                P values of Tmat.
%        DF_Tmat:
%                Degree of freedom of Tmat.
%  TFNBS_f_score:
%                TFNBS-based F values for each edge of interest.
%   TFNBS_f_pval:
%                P values of TFNBS_f_score (FWE-corrected).
%  TFNBS_t_score:
%                Post hoc TFNBS-based T values for each edge of interest.
%   TFNBS_t_pval:
%                P values of TFNBS_t_score (un-corrected, two-tailed).
%     N_sig_edge:
%                The number of edges showing significant group effect.
%   Reg_sig_edge:
%                Nodal indices linked by those edges showing significant
%                group effect.
% TFNBS_f_score_sig_edge:
%                TFNBS-based F values for those edges showing significant
%                group effect.
% TFNBS_f_pval_sig_edge:
%                P values of TFNBS_f_score_sig_edge (FWE-corrected).
% TFNBS_t_score_sig_edge:
%                Post hoc TFNBS-based T values for those edges showing
%                significant group effect.
% TFNBS_t_pval_sig_edge:
%                P values of TFNBS_t_score_sig_edge (un-corrected, two-tailed).
%
% References:
%    Baggio et al., 2018, Statistical inference in brain graphs using threshold-free
%    network-based statistics.
%
% Rui WANG,    IBRR, SCNU, Guangzhou, 2019/10/11, raye_wong@126.com
% Yuping YANG, IBRR, SCNU, Guangzhou, 2019/11/30, yupingyanghvd@gmail.com
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/11/30, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin == 5
    E = 0.5; H = 2;
end

if nargin == 6
    H = 2;
end

Ind_mask = find(triu(Mask_net,1));

N_group = length(Data_path);
N_sub = zeros(N_group,1);
Data = cell(N_group,1);

% reorganize data
for i_group = 1:N_group
    
    cd (Data_path{i_group})
    Mat_file = ls([File_filter{i_group} '*.mat']);
    
    if ~isempty(Mat_file)
        N_sub(i_group) = size(Mat_file,1);
        Data_igroup = zeros(N_sub(i_group),length(Ind_mask));
        
        for i_sub = 1:N_sub(i_group)
            Ind_mat = load(Mat_file(i_sub,:));
            Name_field = fieldnames(Ind_mat);
            Ind_mat = Ind_mat.(Name_field{1});
            Data_igroup(i_sub,:) = Ind_mat(Ind_mask);
        end
    else
        Mat_file = ls([File_filter{i_group} '*.txt']);
        
        N_sub(i_group) = size(Mat_file,1);
        Data_igroup = zeros(N_sub(i_group),length(Ind_mask));
        
        for i_sub = 1:N_sub(i_group)
            Ind_mat = load(Mat_file(i_sub,:));
            Data_igroup(i_sub,:) = Ind_mat(Ind_mask);
        end
    end
    Data{i_group} = Data_igroup;
end

% ancova or anova for real data
if nargin == 8
    N_cova = length(Cova_path);
    
    if N_cova ~= N_group
        error('The number of covariate files is not equal to the number of groups!');
    end
    
    Cova = cell(N_cova,1);
    
    for i_cova = 1:N_cova
        [~,~,ext] = fileparts(Cova_path{i_cova});
        
        switch ext
            case '.txt'
                Cova{i_cova} = load(Cova_path{i_cova});
            case '.mat'
                Ind_cova = load(Cova_path{i_cova});
                Name_field = fieldnames(Ind_cova);
                Cova{i_cova} = Ind_cova.(Name_field{1});
        end
    end
    [Fvec,P_Fvec,DF_Fvec,Post_Tvec, Post_P_Tvec, Post_DF_Tvec] = gretna_ancova1(Data,Cova);
else
    [Fvec,P_Fvec,DF_Fvec,Post_Tvec, Post_P_Tvec, Post_DF_Tvec] = gretna_ancova1(Data);
end

% TFNBS for real data (F test)
Fmat          = zeros(size(Mask_net));
P_Fmat        = zeros(size(Mask_net));
DF_Fmat.group = zeros(size(Mask_net));
DF_Fmat.error = zeros(size(Mask_net));

Fmat(Ind_mask)          = Fvec;
Fmat                    = Fmat + Fmat';
P_Fmat(Ind_mask)        = P_Fvec;
P_Fmat                  = P_Fmat + P_Fmat';
DF_Fmat.group(Ind_mask) = DF_Fvec.group;
DF_Fmat.group           = DF_Fmat.group + DF_Fmat.group';
DF_Fmat.error(Ind_mask) = DF_Fvec.error;
DF_Fmat.error           = DF_Fmat.error + DF_Fmat.error';

F_thr         = linspace(0,max(Fvec),101);
TFNBS_f_score = zeros([size(Mask_net),101]);

for i_thr = 1:101
    F_suprathres                               = Fmat;
    F_suprathres(F_suprathres >= F_thr(i_thr)) = 1;
    F_suprathres(F_suprathres ~= 1)            = 0;
    
    [F_ci, F_n_node] = components(sparse(F_suprathres));
    F_ind_com = find(F_n_node > 1);
    F_n_com = length(F_ind_com);
    
    for i_com = 1:F_n_com
        F_ind_subn = find(F_ci == F_ind_com(i_com));
        
        F_subn = F_suprathres(F_ind_subn,F_ind_subn);
        F_size = sum(sum(F_subn))/2;
        
        F_subn(F_subn~=0) = (F_size^E)*(F_thr(i_thr)^H);
        TFNBS_f_score(F_ind_subn,F_ind_subn,i_thr) = F_subn;
    end
end
TFNBS_f_score = sum(TFNBS_f_score,3);

% TFNBS for real data (post hoc T test)
Tmat    = zeros([size(Mask_net),size(Post_Tvec,1)]);
P_Tmat  = zeros([size(Mask_net),size(Post_Tvec,1)]);
DF_Tmat = zeros([size(Mask_net),size(Post_Tvec,1)]);

for i_post = 1:size(Post_Tvec,1)
    Tmat_i_post    = zeros(size(Mask_net));
    P_Tmat_i_post  = zeros(size(Mask_net));
    DF_Tmat_i_post = zeros(size(Mask_net));
    
    Tmat_i_post(Ind_mask) = Post_Tvec(i_post,:);
    Tmat_i_post           = Tmat_i_post + Tmat_i_post';
    Tmat(:,:,i_post)      = Tmat_i_post;
    
    P_Tmat_i_post(Ind_mask) = Post_P_Tvec(i_post,:);
    P_Tmat_i_post           = P_Tmat_i_post + P_Tmat_i_post';
    P_Tmat(:,:,i_post)      = P_Tmat_i_post;
    
    DF_Tmat_i_post(Ind_mask) = Post_DF_Tvec(i_post,:);
    DF_Tmat_i_post           = DF_Tmat_i_post + DF_Tmat_i_post';
    DF_Tmat(:,:,i_post)      = DF_Tmat_i_post;
end

TFNBS_t_score = zeros([size(Mask_net),size(Post_Tvec,1)]);

for i_post = 1:size(Post_Tvec,1)
    TFNBS_tscore_ipost = zeros([size(Mask_net),101]);
    
    if min(Post_Tvec(i_post,:)) >= 0
        T_thr_pos = linspace(0,max(Post_Tvec(i_post,:)),101);
        
        for i_thr = 1:101
            % positive
            T_suprathres_pos                                       = Tmat(:,:,i_post);
            T_suprathres_pos(T_suprathres_pos >= T_thr_pos(i_thr)) = 1;
            T_suprathres_pos(T_suprathres_pos ~= 1)                = 0;
            
            [T_ci_pos, T_n_node_pos] = components(sparse(T_suprathres_pos));
            
            T_ind_com_pos = find(T_n_node_pos > 1);
            T_n_com_pos   = length(T_ind_com_pos);
            
            for i_com = 1:T_n_com_pos
                T_ind_subn_pos = find(T_ci_pos == T_ind_com_pos(i_com));
                
                T_subn_pos = T_suprathres_pos(T_ind_subn_pos,T_ind_subn_pos);
                T_size_pos = sum(sum(T_subn_pos))/2;
                
                T_subn_pos(T_subn_pos~=0) = (T_size_pos^E)*(T_thr_pos(i_thr)^H);
                TFNBS_tscore_ipost(T_ind_subn_pos,T_ind_subn_pos,i_thr) = T_subn_pos;
            end
        end
        TFNBS_t_score(:,:,i_post) = sum(TFNBS_tscore_ipost,3);
        
    elseif max(Post_Tvec(i_post,:)) <= 0
        T_thr_neg = linspace(min(Post_Tvec(i_post,:)),0,101);
        
        for i_thr = 1:101
            % negative
            T_suprathres_neg                                       = Tmat(:,:,i_post);
            T_suprathres_neg(T_suprathres_neg <= T_thr_neg(i_thr)) = 1;
            T_suprathres_neg(T_suprathres_neg ~= 1)                = 0;
            
            [T_ci_neg, T_n_node_neg] = components(sparse(T_suprathres_neg));
            
            T_ind_com_neg = find(T_n_node_neg > 1);
            T_n_com_neg   = length(T_ind_com_neg);
            
            for i_com = 1:T_n_com_neg
                T_ind_subn_neg = find(T_ci_neg == T_ind_com_neg(i_com));
                
                T_subn_neg = T_suprathres_neg(T_ind_subn_neg,T_ind_subn_neg);
                T_size_neg = sum(sum(T_subn_neg))/2;
                
                T_subn_neg(T_subn_neg~=0) = -(T_size_neg^E)*(T_thr_neg(i_thr)^H);
                TFNBS_tscore_ipost(T_ind_subn_neg,T_ind_subn_neg,i_thr) = T_subn_neg;
            end
        end
        TFNBS_t_score(:,:,i_post) = sum(TFNBS_tscore_ipost,3);
        
    else
        T_thr_pos = linspace(0,max(Post_Tvec(i_post,:)),101);
        T_thr_neg = linspace(min(Post_Tvec(i_post,:)),0,101);
        
        for i_thr = 1:101
            % positive
            T_suprathres_pos                                       = Tmat(:,:,i_post);
            T_suprathres_pos(T_suprathres_pos >= T_thr_pos(i_thr)) = 1;
            T_suprathres_pos(T_suprathres_pos ~= 1)                = 0;
            
            [T_ci_pos, T_n_node_pos] = components(sparse(T_suprathres_pos));
            
            T_ind_com_pos = find(T_n_node_pos > 1);
            T_n_com_pos   = length(T_ind_com_pos);
            
            for i_com = 1:T_n_com_pos
                T_ind_subn_pos = find(T_ci_pos == T_ind_com_pos(i_com));
                
                T_subn_pos = T_suprathres_pos(T_ind_subn_pos,T_ind_subn_pos);
                T_size_pos = sum(sum(T_subn_pos))/2;
                
                T_subn_pos(T_subn_pos~=0) = (T_size_pos^E)*(T_thr_pos(i_thr)^H);
                TFNBS_tscore_ipost(T_ind_subn_pos,T_ind_subn_pos,i_thr) = T_subn_pos;
            end
            
            % negative
            T_suprathres_neg                                       = Tmat(:,:,i_post);
            T_suprathres_neg(T_suprathres_neg <= T_thr_neg(i_thr)) = 1;
            T_suprathres_neg(T_suprathres_neg ~= 1)                = 0;
            
            [T_ci_neg, T_n_node_neg] = components(sparse(T_suprathres_neg));
            
            T_ind_com_neg = find(T_n_node_neg > 1);
            T_n_com_neg   = length(T_ind_com_neg);
            
            for i_com = 1:T_n_com_neg
                T_ind_subn_neg = find(T_ci_neg == T_ind_com_neg(i_com));
                
                T_subn_neg = T_suprathres_neg(T_ind_subn_neg,T_ind_subn_neg);
                T_size_neg = sum(sum(T_subn_neg))/2;
                
                T_subn_neg(T_subn_neg~=0) = -(T_size_neg^E)*(T_thr_neg(i_thr)^H);
                TFNBS_tscore_ipost(T_ind_subn_neg,T_ind_subn_neg,i_thr) = T_subn_neg;
            end
        end
        TFNBS_t_score(:,:,i_post) = sum(TFNBS_tscore_ipost,3);
    end
end

% permutation test
Data = cell2mat(Data);
if nargin == 8
    Cova = cell2mat(Cova);
end

Ind_start = zeros(N_group,1);
Ind_start(1) = 1;
for i_group = 1:N_group-1
    Ind_start(i_group+1) = sum(N_sub(1:i_group))+1;
end

Max_TFNBS_f_score_rand = zeros(M,1);
TFNBS_t_score_rand = zeros([size(Mask_net),size(Post_Tvec,1),M]);

for i_per = 1:M
    Ind_rand = randperm(sum(N_sub));
    Data_rand = cell(N_group, 1);
    
    for i_group = 1:N_group
        Data_rand{i_group} = Data(Ind_rand(Ind_start(i_group):sum(N_sub(1:i_group))),:);
    end
    
    if nargin == 8
        Cova_rand = cell(N_group,1);
        
        for i_group = 1:N_group
            Cova_rand{i_group} = Cova(Ind_rand(Ind_start(i_group):sum(N_sub(1:i_group))),:);
        end
        
        [Fvec_rand,~,~,Post_Tvec_rand,~,~] = gretna_ancova1(Data_rand,Cova_rand);    
    else
        [Fvec_rand,~,~,Post_Tvec_rand,~,~] = gretna_ancova1(Data_rand);
    end
    
    % for F test
    Fmat_rand           = zeros(size(Mask_net));
    Fmat_rand(Ind_mask) = Fvec_rand;
    Fmat_rand           = Fmat_rand + Fmat_rand';
    
    F_thr_rand         = linspace(0, max(Fvec_rand), 101);
    TFNBS_f_score_rand = zeros([size(Mask_net),101]);
    
    for i_thr = 1:101
        F_suprathres_rand                                         = Fmat_rand;
        F_suprathres_rand(F_suprathres_rand >= F_thr_rand(i_thr)) = 1;
        F_suprathres_rand(F_suprathres_rand ~= 1)                 = 0;
        
        [F_ci_rand, F_n_node_rand] = components(sparse(F_suprathres_rand));
        
        F_ind_com_rand = find(F_n_node_rand > 1);
        F_n_com_rand = length(F_ind_com_rand);
        
        for i_com = 1:F_n_com_rand
            F_ind_subn_rand = find(F_ci_rand == F_ind_com_rand(i_com));
            
            F_subn_rand = F_suprathres_rand(F_ind_subn_rand,F_ind_subn_rand);
            F_size_rand = sum(sum(F_subn_rand))/2;
            
            F_subn_rand(F_subn_rand~=0) = (F_size_rand^E)*(F_thr_rand(i_thr)^H);
            TFNBS_f_score_rand(F_ind_subn_rand,F_ind_subn_rand,i_thr) = F_subn_rand;
        end
    end
    TFNBS_f_score_rand = sum(TFNBS_f_score_rand,3);
    Max_TFNBS_f_score_rand(i_per) = max(TFNBS_f_score_rand(:));
    
    % for post hoc T test
    for i_post = 1:size(Post_Tvec_rand,1)
        Tmat_rand           = zeros(size(Mask_net));
        Tmat_rand(Ind_mask) = Post_Tvec_rand(i_post,:);
        Tmat_rand           = Tmat_rand + Tmat_rand';
        
        TFNBS_t_score_ipost_rand = zeros([size(Mask_net),101]);
        
        if min(Post_Tvec_rand(i_post,:)) >= 0
            T_thr_pos_rand = linspace(0,max(Post_Tvec_rand(i_post,:)),101);
            
            for i_thr = 1:101
                % positive
                T_suprathres_pos_rand                                                 = Tmat_rand;
                T_suprathres_pos_rand(T_suprathres_pos_rand >= T_thr_pos_rand(i_thr)) = 1;
                T_suprathres_pos_rand(T_suprathres_pos_rand ~= 1)                     = 0;
                
                [T_ci_pos_rand, T_n_node_pos_rand] = components(sparse(T_suprathres_pos_rand));
                
                T_ind_com_pos_rand = find(T_n_node_pos_rand > 1);
                T_n_com_pos_rand   = length(T_ind_com_pos_rand);
                
                for i_com = 1:T_n_com_pos_rand
                    T_ind_subn_pos_rand = find(T_ci_pos_rand == T_ind_com_pos_rand(i_com));
                    
                    T_subn_pos_rand = T_suprathres_pos_rand(T_ind_subn_pos_rand,T_ind_subn_pos_rand);
                    T_size_pos_rand = sum(sum(T_subn_pos_rand))/2;
                    
                    T_subn_pos_rand(T_subn_pos_rand~=0) = (T_size_pos_rand^E)*(T_thr_pos_rand(i_thr)^H);
                    TFNBS_t_score_ipost_rand(T_ind_subn_pos_rand,T_ind_subn_pos_rand,i_thr) = T_subn_pos_rand;
                end
            end
            TFNBS_t_score_rand(:,:,i_post,i_per) = sum(TFNBS_t_score_ipost_rand,3);
            
        elseif max(Post_Tvec_rand(i_post,:)) <= 0
            T_thr_neg_rand = linspace(min(Post_Tvec_rand(i_post,:)),0,101);
            
            for i_thr = 1:101
                % negative
                T_suprathres_neg_rand                                                 = Tmat_rand;
                T_suprathres_neg_rand(T_suprathres_neg_rand <= T_thr_neg_rand(i_thr)) = 1;
                T_suprathres_neg_rand(T_suprathres_neg_rand ~= 1)                     = 0;
                
                [T_ci_neg_rand, T_n_node_neg_rand] = components(sparse(T_suprathres_neg_rand));
                
                T_ind_com_neg_rand = find(T_n_node_neg_rand > 1);
                T_n_com_neg_rand   = length(T_ind_com_neg_rand);
                
                for i_com = 1:T_n_com_neg_rand
                    T_ind_subn_neg_rand = find(T_ci_neg_rand == T_ind_com_neg_rand(i_com));
                    
                    T_subn_neg_rand = T_suprathres_neg_rand(T_ind_subn_neg_rand,T_ind_subn_neg_rand);
                    T_size_neg_rand = sum(sum(T_subn_neg_rand))/2;
                    
                    T_subn_neg_rand(T_subn_neg_rand~=0) = -(T_size_neg_rand^E)*(T_thr_neg_rand(i_thr)^H);
                    TFNBS_t_score_ipost_rand(T_ind_subn_neg_rand,T_ind_subn_neg_rand,i_thr) = T_subn_neg_rand;
                end
            end
            TFNBS_t_score_rand(:,:,i_post,i_per) = sum(TFNBS_t_score_ipost_rand,3);
            
        else
            T_thr_pos_rand = linspace(0,max(Post_Tvec_rand(i_post,:)),101);
            T_thr_neg_rand = linspace(min(Post_Tvec_rand(i_post,:)),0,101);
            
            for i_thr = 1:101
                % positive
                T_suprathres_pos_rand                                                 = Tmat_rand;
                T_suprathres_pos_rand(T_suprathres_pos_rand >= T_thr_pos_rand(i_thr)) = 1;
                T_suprathres_pos_rand(T_suprathres_pos_rand ~= 1)                     = 0;
                
                [T_ci_pos_rand, T_n_node_pos_rand] = components(sparse(T_suprathres_pos_rand));
                
                T_ind_com_pos_rand = find(T_n_node_pos_rand > 1);
                T_n_com_pos_rand   = length(T_ind_com_pos_rand);
                
                for i_com = 1:T_n_com_pos_rand
                    T_ind_subn_pos_rand = find(T_ci_pos_rand == T_ind_com_pos_rand(i_com));
                    
                    T_subn_pos_rand = T_suprathres_pos_rand(T_ind_subn_pos_rand,T_ind_subn_pos_rand);
                    T_size_pos_rand = sum(sum(T_subn_pos_rand))/2;
                    
                    T_subn_pos_rand(T_subn_pos_rand~=0) = (T_size_pos_rand^E)*(T_thr_pos_rand(i_thr)^H);
                    TFNBS_t_score_ipost_rand(T_ind_subn_pos_rand,T_ind_subn_pos_rand,i_thr) = T_subn_pos_rand;
                end
                
                % negative
                T_suprathres_neg_rand                                                 = Tmat_rand;
                T_suprathres_neg_rand(T_suprathres_neg_rand <= T_thr_neg_rand(i_thr)) = 1;
                T_suprathres_neg_rand(T_suprathres_neg_rand ~= 1)                     = 0;
                
                [T_ci_neg_rand, T_n_node_neg_rand] = components(sparse(T_suprathres_neg_rand));
                
                T_ind_com_neg_rand = find(T_n_node_neg_rand > 1);
                T_n_com_neg_rand   = length(T_ind_com_neg_rand);
                
                for i_com = 1:T_n_com_neg_rand
                    T_ind_subn_neg_rand = find(T_ci_neg_rand == T_ind_com_neg_rand(i_com));
                    
                    T_subn_neg_rand = T_suprathres_neg_rand(T_ind_subn_neg_rand,T_ind_subn_neg_rand);
                    T_size_neg_rand = sum(sum(T_subn_neg_rand))/2;
                    
                    T_subn_neg_rand(T_subn_neg_rand~=0) = -(T_size_neg_rand^E)*(T_thr_neg_rand(i_thr)^H);
                    TFNBS_t_score_ipost_rand(T_ind_subn_neg_rand,T_ind_subn_neg_rand,i_thr) = T_subn_neg_rand;
                end
            end
            TFNBS_t_score_rand(:,:,i_post,i_per) = sum(TFNBS_t_score_ipost_rand,3);
        end
    end
    
    fprintf('%d permutation test(s) have been done\n', i_per);
    
end

% calculate FWE-corrected p value (TFNBS_fscore)
TFNBS_f_pval = zeros(size(Mask_net));
for i_edge = 1:length(Ind_mask)
    TFNBS_f_pval(Ind_mask(i_edge)) = sum(Max_TFNBS_f_score_rand > TFNBS_f_score(Ind_mask(i_edge)))/(M+1);
end
TFNBS_f_pval = TFNBS_f_pval + TFNBS_f_pval';

% calculate uncorrected, two-tailed p value (post hoc TFNBS_tscore)
TFNBS_t_pval = zeros([size(Mask_net),size(Post_Tvec,1)]);
for i_post = 1:size(Post_Tvec,1)
    TFNBS_t_pval_i_post = zeros(size(Mask_net));
    for i_edge = 1:length(Ind_mask)
        [Row,Col] = ind2sub(size(Mask_net),Ind_mask(i_edge));
        TFNBS_t_score_rand_i_post_i_edge = squeeze(TFNBS_t_score_rand(Row,Col,i_post,:));
        TFNBS_t_pval_i_post(Row,Col) = sum(abs(TFNBS_t_score_rand_i_post_i_edge)>abs(TFNBS_t_score(Row,Col,i_post)))/(M+1);
    end
    TFNBS_t_pval_i_post = TFNBS_t_pval_i_post + TFNBS_t_pval_i_post';
    TFNBS_t_pval(:,:,i_post) = TFNBS_t_pval_i_post;
end

% determine edges showing significant group effect
Ind_sig_edge = find(TFNBS_f_pval(Ind_mask) < Pthr);
N_sig_edge = length(Ind_sig_edge);
[Reg_sig_edge(:,1),Reg_sig_edge(:,2)] = ind2sub(size(Mask_net),Ind_mask(Ind_sig_edge));

TFNBS_f_score_sig_edge = TFNBS_f_score(Ind_mask(Ind_sig_edge));
TFNBS_f_pval_sig_edge  = TFNBS_f_pval(Ind_mask(Ind_sig_edge));

TFNBS_t_score_sig_edge = zeros(N_sig_edge,size(Post_Tvec,1));
TFNBS_t_pval_sig_edge  = zeros(N_sig_edge,size(Post_Tvec,1));
for i_post = 1:size(Post_Tvec,1)
    TFNBS_t_score_i_post = TFNBS_t_score(:,:,i_post);
    TFNBS_t_pval_i_post  = TFNBS_t_pval(:,:,i_post);
    for i_edge = 1:N_sig_edge
        TFNBS_t_score_sig_edge(:,i_post) = TFNBS_t_score_i_post(Ind_mask(Ind_sig_edge));
        TFNBS_t_pval_sig_edge(:,i_post)  = TFNBS_t_pval_i_post(Ind_mask(Ind_sig_edge));
    end
end

TFNBS_Result.Fmat    = Fmat;
TFNBS_Result.P_Fmat  = P_Fmat;
TFNBS_Result.DF_Fmat = DF_Fmat;
TFNBS_Result.Tmat    = Tmat;
TFNBS_Result.P_Tmat  = P_Tmat;
TFNBS_Result.DF_Tmat = DF_Tmat;

TFNBS_Result.TFNBS_f_score = TFNBS_f_score;
TFNBS_Result.TFNBS_f_pval  = TFNBS_f_pval;
TFNBS_Result.TFNBS_t_score = TFNBS_t_score;
TFNBS_Result.TFNBS_t_pval  = TFNBS_t_pval;

TFNBS_Result.N_sig_edge             = N_sig_edge;
TFNBS_Result.Reg_sig_edge           = Reg_sig_edge;
TFNBS_Result.TFNBS_f_score_sig_edge = TFNBS_f_score_sig_edge;
TFNBS_Result.TFNBS_f_pval_sig_edge  = TFNBS_f_pval_sig_edge;
TFNBS_Result.TFNBS_t_score_sig_edge = TFNBS_t_score_sig_edge;
TFNBS_Result.TFNBS_t_pval_sig_edge  = TFNBS_t_pval_sig_edge;

return