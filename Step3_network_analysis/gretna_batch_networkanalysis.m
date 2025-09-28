function [Net_para, Node_para] = gretna_batch_networkanalysis(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype)

%==========================================================================
% This function is used to calculate several global and nodal network
% metrics for both binary and weighted networks over a continuous threshold
% range for multiple subjects.
%
%
% Syntax: function [Net_para, Node_para] = gretna_batch_networkanalysis(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype)
%
% Inputs:
%          Mat_path:
%                   The directory where individual matrices are sorted.
%       File_filter:
%                   The prefix of those matrices.
%              Thr1:
%                   The lower bound of the threshold range.
%              Thr2:
%                   The upper bound of the threhsold range.
%             Delta:
%                   The interval of the threshold range.
%            N_rand:
%                   The number of random networks.
%             Stype:
%                   'pos' for only positive elements;
%                   'neg' for only negative elements using absolute values;
%                   'abs' for all elements using absolute values.
%             Ttype:
%                   'r' for correlation threshold;
%                   's' for sparsity threshold;
%                   'k' for edge threshold.
%
% Outputs:
%          Net_para:
%                   Global network metrics.
%         Node_para:
%                   Nodal network metrics.
%
%               .bin:
%                   Binary networks.
%               .wei:
%                   Weighted networks.
%
%               .cp:
%                   Clustering coefficient.
%               .lp:
%                   Shortest path length.
%             .loce:
%                   Local efficiency.
%               .ge:
%                   Global efficiency.
%              .deg:
%                   Degree.
%               .bw:
%                   Betweennesses.
%               .bg:
%                   Bridging centrality.
%              .clo:
%                   Closeness centrality.
%              .com:
%                   Communicability centrality.
%              .eig:
%                   Eigenvector centrality.
%             .hind:
%                   H-index centrality.
%              .inf:
%                   Information centrality.
%             .katz:
%                   Katz centrality.
%              .lap:
%                   Laplacian centrality.
%              .lev:
%                   Leverage centrality.
%             .page:
%                   Pagerank centrality.
%              .ass:
%                   Assortativity.
%             .beta:
%                   Hierarchy.
%                .q:
%                   Modularity.
%
%          .XXXrand:
%                   Each network metric of random networks.
%         .XXXratio:
%                   Ratio of each network metric between real and random
%                   networks.
%              aXXX:
%                   Area under curve of each network metric.
%
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/10/31, jinhui.wang.1982@gmail.com
%==========================================================================

Thres = Thr1:Delta:Thr2;

cd(Mat_path)
Mats = ls([File_filter '*.mat']);

for isub = 1:size(Mats,1) % subjects
    
    fprintf('------------------------------------------------------------------------------------\n')
    fprintf('Calculating network parameters for %s\n', Mats(isub,:));
    fprintf('------------------------------------------------------------------------------------\n')
    
    Matrix = load(deblank(Mats(isub,:)));
    Field_name = fieldnames(Matrix);
    Matrix = Matrix.(Field_name{1});
    
    for ithres = 1:length(Thres) % thresholds
        
        [Bin, Wei] = gretna_R2b_MST(Matrix,Stype,Ttype,Thres(ithres));
        
        % clustering coefficient
        [Net_para.bin.cp(ithres,isub), Node_para.bin.cp(ithres,isub,:)] = gretna_node_clustcoeff(Bin);
        [Net_para.wei.cp(ithres,isub), Node_para.wei.cp(ithres,isub,:)] = gretna_node_clustcoeff_weight(Wei,1);
        
        % shortest path length
        [Net_para.bin.lp(ithres,isub), Node_para.bin.lp(ithres,isub,:)] = gretna_node_shortestpathlength(Bin);
        [Net_para.wei.lp(ithres,isub), Node_para.wei.lp(ithres,isub,:)] = gretna_node_shortestpathlength_weight(Wei);
        
        % local efficiency
        [Net_para.bin.loce(ithres,isub), Node_para.bin.loce(ithres,isub,:)] = gretna_node_local_efficiency(Bin);
        [Net_para.wei.loce(ithres,isub), Node_para.wei.loce(ithres,isub,:)] = gretna_node_local_efficiency_weight(Wei);
        
        % global efficiency
        [Net_para.bin.ge(ithres,isub), Node_para.bin.ge(ithres,isub,:)] = gretna_node_global_efficiency(Bin);
        [Net_para.wei.ge(ithres,isub), Node_para.wei.ge(ithres,isub,:)] = gretna_node_global_efficiency_weight(Wei);
        
        % degree
        [Net_para.bin.deg(ithres,isub), Node_para.bin.deg(ithres,isub,:)] = gretna_node_degree(Bin);
        [Net_para.wei.deg(ithres,isub), Node_para.wei.deg(ithres,isub,:)] = gretna_node_degree_weight(Wei);
        
        % betweenness
        [Net_para.bin.bw(ithres,isub), Node_para.bin.bw(ithres,isub,:)] = gretna_node_betweenness(Bin);
        [Net_para.wei.bw(ithres,isub), Node_para.wei.bw(ithres,isub,:)] = gretna_node_betweenness_weight(Wei);
        
        % bridging centrality
        [Net_para.bin.bg(ithres,isub), Node_para.bin.bg(ithres,isub,:)] = gretna_node_bridging_centrality(Bin);
        [Net_para.wei.bg(ithres,isub), Node_para.wei.bg(ithres,isub,:)] = gretna_node_bridging_centrality_weight(Wei);
        
        % closeness centrality
        [Net_para.bin.clo(ithres,isub), Node_para.bin.clo(ithres,isub,:)] = gretna_node_closeness_centrality(Bin);
        [Net_para.wei.clo(ithres,isub), Node_para.wei.clo(ithres,isub,:)] = gretna_node_closeness_centrality_weight(Wei);
        
        % communicability centrality
        [Net_para.bin.com(ithres,isub), Node_para.bin.com(ithres,isub,:)] = gretna_node_communicability_betweenness(Bin);
        [Net_para.wei.com(ithres,isub), Node_para.wei.com(ithres,isub,:)] = gretna_node_communicability_betweenness_weight(Wei);
        
        % eigenvector centrality
        [Net_para.bin.eig(ithres,isub), Node_para.bin.eig(ithres,isub,:)] = gretna_node_eigenvector(Bin);
        [Net_para.wei.eig(ithres,isub), Node_para.wei.eig(ithres,isub,:)] = gretna_node_eigenvector_weight(Wei);
        
        % h-index centrality
        [Net_para.bin.hind(ithres,isub), Node_para.bin.hind(ithres,isub,:)] = gretna_node_hindex(Bin);
        [Net_para.wei.hind(ithres,isub), Node_para.wei.hind(ithres,isub,:)] = gretna_node_hindex_weight(Wei);
        
        % information centrality
        [Net_para.bin.inf(ithres,isub), Node_para.bin.inf(ithres,isub,:)] = gretna_node_information_centrality(Bin);
        [Net_para.wei.inf(ithres,isub), Node_para.wei.inf(ithres,isub,:)] = gretna_node_information_centrality_weight(Wei);
        
        % katz centrality
        [Net_para.bin.katz(ithres,isub), Node_para.bin.katz(ithres,isub,:)] = gretna_node_katz_centrality(Bin);
        [Net_para.wei.katz(ithres,isub), Node_para.wei.katz(ithres,isub,:)] = gretna_node_katz_centrality_weight(Wei);
        
        % laplacian centrality
        [Net_para.bin.lap(ithres,isub), Node_para.bin.lap(ithres,isub,:)] = gretna_node_laplacian_centrality(Bin);
        [Net_para.wei.lap(ithres,isub), Node_para.wei.lap(ithres,isub,:)] = gretna_node_laplacian_centrality_weight(Wei);
        
        % leverage centrality
        [Net_para.bin.lev(ithres,isub), Node_para.bin.lev(ithres,isub,:)] = gretna_node_leverage_centrality(Bin);
        [Net_para.wei.lev(ithres,isub), Node_para.wei.lev(ithres,isub,:)] = gretna_node_leverage_centrality_weight(Wei);
        
        % pagerank centrality
        [Net_para.bin.page(ithres,isub), Node_para.bin.page(ithres,isub,:)] = gretna_node_pagerank(Bin);
        [Net_para.wei.page(ithres,isub), Node_para.wei.page(ithres,isub,:)] = gretna_node_pagerank_weight(Wei);
        
        % assortativity
        [Net_para.bin.ass(ithres,isub)] = gretna_assortativity(Bin);
        [Net_para.wei.ass(ithres,isub)] = gretna_assortativity_weight(Wei);
        
        % hierarchy
        [Net_para.bin.beta(ithres,isub)] = gretna_hierarchy(Bin);
        [Net_para.wei.beta(ithres,isub)] = gretna_hierarchy_weight(Wei);
        
        % modularity
        [~,Net_para.bin.q(ithres,isub)] = modularity_und(Bin);
        [~,Net_para.wei.q(ithres,isub)] = modularity_und(Wei);
        
        for irand = 1:N_rand % random networks
            
            Bin_rand = gretna_gen_random_network1(Bin);
            Wei_rand = gretna_gen_random_network1_weight(Wei);
            
            % clustering coefficient
            [Net_para.bin.cprand(ithres,isub,irand), Node_para.bin.cprand(ithres,isub,:,irand)] = gretna_node_clustcoeff(Bin_rand);
            [Net_para.wei.cprand(ithres,isub,irand), Node_para.wei.cprand(ithres,isub,:,irand)] = gretna_node_clustcoeff_weight(Wei_rand,1);
            
            % shortest path length
            [Net_para.bin.lprand(ithres,isub,irand), Node_para.bin.lprand(ithres,isub,:,irand)] = gretna_node_shortestpathlength(Bin_rand);
            [Net_para.wei.lprand(ithres,isub,irand), Node_para.wei.lprand(ithres,isub,:,irand)] = gretna_node_shortestpathlength_weight(Wei_rand);
            
            % local efficiency
            [Net_para.bin.locerand(ithres,isub,irand), Node_para.bin.locerand(ithres,isub,:,irand)] = gretna_node_local_efficiency(Bin_rand);
            [Net_para.wei.locerand(ithres,isub,irand), Node_para.wei.locerand(ithres,isub,:,irand)] = gretna_node_local_efficiency_weight(Wei_rand);
            
            % global efficiency
            [Net_para.bin.gerand(ithres,isub,irand), Node_para.bin.gerand(ithres,isub,:,irand)] = gretna_node_global_efficiency(Bin_rand);
            [Net_para.wei.gerand(ithres,isub,irand), Node_para.wei.gerand(ithres,isub,:,irand)] = gretna_node_global_efficiency_weight(Wei_rand);
            
            % degree
            [Net_para.bin.degrand(ithres,isub,irand), Node_para.bin.degrand(ithres,isub,:,irand)] = gretna_node_degree(Bin_rand);
            [Net_para.wei.degrand(ithres,isub,irand), Node_para.wei.degrand(ithres,isub,:,irand)] = gretna_node_degree_weight(Wei_rand);
            
            % betweenness
            [Net_para.bin.bwrand(ithres,isub,irand), Node_para.bin.bwrand(ithres,isub,:,irand)]  = gretna_node_betweenness(Bin_rand);
            [Net_para.wei.bwrand(ithres,isub,irand), Node_para.wei.bwrand(ithres,isub,:,irand)]  = gretna_node_betweenness_weight(Wei_rand);
            
            % bridging centrality
            [Net_para.bin.bgrand(ithres,isub,irand), Node_para.bin.bgrand(ithres,isub,:,irand)] = gretna_node_bridging_centrality(Bin_rand);
            [Net_para.wei.bgrand(ithres,isub,irand), Node_para.wei.bgrand(ithres,isub,:,irand)] = gretna_node_bridging_centrality_weight(Wei_rand);
            
            % closeness centrality
            [Net_para.bin.clorand(ithres,isub,irand), Node_para.bin.clorand(ithres,isub,:,irand)] = gretna_node_closeness_centrality(Bin_rand);
            [Net_para.wei.clorand(ithres,isub,irand), Node_para.wei.clorand(ithres,isub,:,irand)] = gretna_node_closeness_centrality_weight(Wei_rand);
            
            % communicability centrality
            [Net_para.bin.comrand(ithres,isub,irand), Node_para.bin.comrand(ithres,isub,:,irand)] = gretna_node_communicability_betweenness(Bin_rand);
            [Net_para.wei.comrand(ithres,isub,irand), Node_para.wei.comrand(ithres,isub,:,irand)] = gretna_node_communicability_betweenness_weight(Wei_rand);
            
            % eigenvector centrality
            [Net_para.bin.eigrand(ithres,isub,irand), Node_para.bin.eigrand(ithres,isub,:,irand)] = gretna_node_eigenvector(Bin_rand);
            [Net_para.wei.eigrand(ithres,isub,irand), Node_para.wei.eigrand(ithres,isub,:,irand)] = gretna_node_eigenvector_weight(Wei_rand);
            
            % h-index centrality
            [Net_para.bin.hindrand(ithres,isub,irand), Node_para.bin.hindrand(ithres,isub,:,irand)] = gretna_node_hindex(Bin_rand);
            [Net_para.wei.hindrand(ithres,isub,irand), Node_para.wei.hindrand(ithres,isub,:,irand)] = gretna_node_hindex_weight(Wei_rand);
            
            % information centrality
            [Net_para.bin.infrand(ithres,isub,irand), Node_para.bin.infrand(ithres,isub,:,irand)] = gretna_node_information_centrality(Bin_rand);
            [Net_para.wei.infrand(ithres,isub,irand), Node_para.wei.infrand(ithres,isub,:,irand)] = gretna_node_information_centrality_weight(Wei_rand);
            
            % katz centrality
            [Net_para.bin.katzrand(ithres,isub,irand), Node_para.bin.katzrand(ithres,isub,:,irand)] = gretna_node_katz_centrality(Bin_rand);
            [Net_para.wei.katzrand(ithres,isub,irand), Node_para.wei.katzrand(ithres,isub,:,irand)] = gretna_node_katz_centrality_weight(Wei_rand);
            
            % laplacian centrality
            [Net_para.bin.laprand(ithres,isub,irand), Node_para.bin.laprand(ithres,isub,:,irand)] = gretna_node_laplacian_centrality(Bin_rand);
            [Net_para.wei.laprand(ithres,isub,irand), Node_para.wei.laprand(ithres,isub,:,irand)] = gretna_node_laplacian_centrality_weight(Wei_rand);
            
            % leverage centrality
            [Net_para.bin.levrand(ithres,isub,irand), Node_para.bin.levrand(ithres,isub,:,irand)] = gretna_node_leverage_centrality(Bin_rand);
            [Net_para.wei.levrand(ithres,isub,irand), Node_para.wei.levrand(ithres,isub,:,irand)] = gretna_node_leverage_centrality_weight(Wei_rand);
            
            % pagerank centrality
            [Net_para.bin.pagerand(ithres,isub,irand), Node_para.bin.pagerand(ithres,isub,:,irand)] = gretna_node_pagerank(Bin_rand);
            [Net_para.wei.pagerand(ithres,isub,irand), Node_para.wei.pagerand(ithres,isub,:,irand)] = gretna_node_pagerank_weight(Wei_rand);
            
            % assortativity
            [Net_para.bin.assrand(ithres,isub,irand)] = gretna_assortativity(Bin_rand);
            [Net_para.wei.assrand(ithres,isub,irand)] = gretna_assortativity_weight(Wei_rand);
            
            % hierarchy
            [Net_para.bin.betarand(ithres,isub,irand)] = gretna_hierarchy(Bin_rand);
            [Net_para.wei.betarand(ithres,isub,irand)] = gretna_hierarchy_weight(Wei_rand);
            
            % modularity
            [~,Net_para.bin.qrand(ithres,isub,irand)] = modularity_und(Bin_rand);
            [~,Net_para.wei.qrand(ithres,isub,irand)] = modularity_und(Wei_rand);
            
        end % random networks
        
        fprintf ('Threshold = %.3f ...... done\n', Thres(ithres));
        
    end % thresholds
    
end % subjects

% normalization
Net_para.bin.cpratio   = Net_para.bin.cp   ./ nanmean(Net_para.bin.cprand,3);
Net_para.bin.lpratio   = Net_para.bin.lp   ./ nanmean(Net_para.bin.lprand,3);
Net_para.bin.loceratio = Net_para.bin.loce ./ nanmean(Net_para.bin.locerand,3);
Net_para.bin.geratio   = Net_para.bin.ge   ./ nanmean(Net_para.bin.gerand,3);
Net_para.bin.degratio  = Net_para.bin.deg  ./ nanmean(Net_para.bin.degrand,3);
Net_para.bin.bwratio   = Net_para.bin.bw   ./ nanmean(Net_para.bin.bwrand,3);
Net_para.bin.bgratio   = Net_para.bin.bg   ./ nanmean(Net_para.bin.bgrand,3);
Net_para.bin.cloratio  = Net_para.bin.clo  ./ nanmean(Net_para.bin.clorand,3);
Net_para.bin.comratio  = Net_para.bin.com  ./ nanmean(Net_para.bin.comrand,3);
Net_para.bin.eigratio  = Net_para.bin.eig  ./ nanmean(Net_para.bin.eigrand,3);
Net_para.bin.hindratio = Net_para.bin.hind ./ nanmean(Net_para.bin.hindrand,3);
Net_para.bin.infratio  = Net_para.bin.inf  ./ nanmean(Net_para.bin.infrand,3);
Net_para.bin.katzratio = Net_para.bin.katz ./ nanmean(Net_para.bin.katzrand,3);
Net_para.bin.lapratio  = Net_para.bin.lap  ./ nanmean(Net_para.bin.laprand,3);
Net_para.bin.levratio  = Net_para.bin.lev  ./ nanmean(Net_para.bin.levrand,3);
Net_para.bin.pageratio = Net_para.bin.page ./ nanmean(Net_para.bin.pagerand,3);
Net_para.bin.assratio  = Net_para.bin.ass  ./ nanmean(Net_para.bin.assrand,3);
Net_para.bin.betaratio = Net_para.bin.beta ./ nanmean(Net_para.bin.betarand,3);
Net_para.bin.qratio    = Net_para.bin.q    ./ nanmean(Net_para.bin.qrand,3);

Node_para.bin.cpratio   = Node_para.bin.cp   ./ nanmean(Node_para.bin.cprand,4);
Node_para.bin.lpratio   = Node_para.bin.lp   ./ nanmean(Node_para.bin.lprand,4);
Node_para.bin.loceratio = Node_para.bin.loce ./ nanmean(Node_para.bin.locerand,4);
Node_para.bin.geratio   = Node_para.bin.ge   ./ nanmean(Node_para.bin.gerand,4);
Node_para.bin.degratio  = Node_para.bin.deg  ./ nanmean(Node_para.bin.degrand,4);
Node_para.bin.bwratio   = Node_para.bin.bw   ./ nanmean(Node_para.bin.bwrand,4);
Node_para.bin.bgratio   = Node_para.bin.bg   ./ nanmean(Node_para.bin.bgrand,4);
Node_para.bin.cloratio  = Node_para.bin.clo  ./ nanmean(Node_para.bin.clorand,4);
Node_para.bin.comratio  = Node_para.bin.com  ./ nanmean(Node_para.bin.comrand,4);
Node_para.bin.eigratio  = Node_para.bin.eig  ./ nanmean(Node_para.bin.eigrand,4);
Node_para.bin.hindratio = Node_para.bin.hind ./ nanmean(Node_para.bin.hindrand,4);
Node_para.bin.infratio  = Node_para.bin.inf  ./ nanmean(Node_para.bin.infrand,4);
Node_para.bin.katzratio = Node_para.bin.katz ./ nanmean(Node_para.bin.katzrand,4);
Node_para.bin.lapratio  = Node_para.bin.lap  ./ nanmean(Node_para.bin.laprand,4);
Node_para.bin.levratio  = Node_para.bin.lev  ./ nanmean(Node_para.bin.levrand,4);
Node_para.bin.pageratio = Node_para.bin.page ./ nanmean(Node_para.bin.pagerand,4);


Net_para.wei.cpratio   = Net_para.wei.cp   ./ nanmean(Net_para.wei.cprand,3);
Net_para.wei.lpratio   = Net_para.wei.lp   ./ nanmean(Net_para.wei.lprand,3);
Net_para.wei.loceratio = Net_para.wei.loce ./ nanmean(Net_para.wei.locerand,3);
Net_para.wei.geratio   = Net_para.wei.ge   ./ nanmean(Net_para.wei.gerand,3);
Net_para.wei.degratio  = Net_para.wei.deg  ./ nanmean(Net_para.wei.degrand,3);
Net_para.wei.bwratio   = Net_para.wei.bw   ./ nanmean(Net_para.wei.bwrand,3);
Net_para.wei.bgratio   = Net_para.wei.bg   ./ nanmean(Net_para.wei.bgrand,3);
Net_para.wei.cloratio  = Net_para.wei.clo  ./ nanmean(Net_para.wei.clorand,3);
Net_para.wei.comratio  = Net_para.wei.com  ./ nanmean(Net_para.wei.comrand,3);
Net_para.wei.eigratio  = Net_para.wei.eig  ./ nanmean(Net_para.wei.eigrand,3);
Net_para.wei.hindratio = Net_para.wei.hind ./ nanmean(Net_para.wei.hindrand,3);
Net_para.wei.infratio  = Net_para.wei.inf  ./ nanmean(Net_para.wei.infrand,3);
Net_para.wei.katzratio = Net_para.wei.katz ./ nanmean(Net_para.wei.katzrand,3);
Net_para.wei.lapratio  = Net_para.wei.lap  ./ nanmean(Net_para.wei.laprand,3);
Net_para.wei.levratio  = Net_para.wei.lev  ./ nanmean(Net_para.wei.levrand,3);
Net_para.wei.pageratio = Net_para.wei.page ./ nanmean(Net_para.wei.pagerand,3);
Net_para.wei.assratio  = Net_para.wei.ass  ./ nanmean(Net_para.wei.assrand,3);
Net_para.wei.betaratio = Net_para.wei.beta ./ nanmean(Net_para.wei.betarand,3);
Net_para.wei.qratio    = Net_para.wei.q    ./ nanmean(Net_para.wei.qrand,3);

Node_para.wei.cpratio   = Node_para.wei.cp   ./ nanmean(Node_para.wei.cprand,4);
Node_para.wei.lpratio   = Node_para.wei.lp   ./ nanmean(Node_para.wei.lprand,4);
Node_para.wei.loceratio = Node_para.wei.loce ./ nanmean(Node_para.wei.locerand,4);
Node_para.wei.geratio   = Node_para.wei.ge   ./ nanmean(Node_para.wei.gerand,4);
Node_para.wei.degratio  = Node_para.wei.deg  ./ nanmean(Node_para.wei.degrand,4);
Node_para.wei.bwratio   = Node_para.wei.bw   ./ nanmean(Node_para.wei.bwrand,4);
Node_para.wei.bgratio   = Node_para.wei.bg   ./ nanmean(Node_para.wei.bgrand,4);
Node_para.wei.cloratio  = Node_para.wei.clo  ./ nanmean(Node_para.wei.clorand,4);
Node_para.wei.comratio  = Node_para.wei.com  ./ nanmean(Node_para.wei.comrand,4);
Node_para.wei.eigratio  = Node_para.wei.eig  ./ nanmean(Node_para.wei.eigrand,4);
Node_para.wei.hindratio = Node_para.wei.hind ./ nanmean(Node_para.wei.hindrand,4);
Node_para.wei.infratio  = Node_para.wei.inf  ./ nanmean(Node_para.wei.infrand,4);
Node_para.wei.katzratio = Node_para.wei.katz ./ nanmean(Node_para.wei.katzrand,4);
Node_para.wei.lapratio  = Node_para.wei.lap  ./ nanmean(Node_para.wei.laprand,4);
Node_para.wei.levratio  = Node_para.wei.lev  ./ nanmean(Node_para.wei.levrand,4);
Node_para.wei.pageratio = Node_para.wei.page ./ nanmean(Node_para.wei.pagerand,4);

% area under curve
for isub = 1:size(Mats,1)
    Net_para.bin.acp(isub,1)   = gretna_auc(Net_para.bin.cp(:,isub),Delta);
    Net_para.bin.alp(isub,1)   = gretna_auc(Net_para.bin.lp(:,isub),Delta);
    Net_para.bin.aloce(isub,1) = gretna_auc(Net_para.bin.loce(:,isub),Delta);
    Net_para.bin.age(isub,1)   = gretna_auc(Net_para.bin.ge(:,isub),Delta);
    Net_para.bin.adeg(isub,1)  = gretna_auc(Net_para.bin.deg(:,isub),Delta);
    Net_para.bin.abw(isub,1)   = gretna_auc(Net_para.bin.bw(:,isub),Delta);
    Net_para.bin.abg(isub,1)   = gretna_auc(Net_para.bin.bg(:,isub),Delta);
    Net_para.bin.aclo(isub,1)  = gretna_auc(Net_para.bin.clo(:,isub),Delta);
    Net_para.bin.acom(isub,1)  = gretna_auc(Net_para.bin.com(:,isub),Delta);
    Net_para.bin.aeig(isub,1)  = gretna_auc(Net_para.bin.eig(:,isub),Delta);
    Net_para.bin.ahind(isub,1) = gretna_auc(Net_para.bin.hind(:,isub),Delta);
    Net_para.bin.ainf(isub,1)  = gretna_auc(Net_para.bin.inf(:,isub),Delta);
    Net_para.bin.akatz(isub,1) = gretna_auc(Net_para.bin.katz(:,isub),Delta);
    Net_para.bin.alap(isub,1)  = gretna_auc(Net_para.bin.lap(:,isub),Delta);
    Net_para.bin.alev(isub,1)  = gretna_auc(Net_para.bin.lev(:,isub),Delta);
    Net_para.bin.apage(isub,1) = gretna_auc(Net_para.bin.page(:,isub),Delta);
    Net_para.bin.aass(isub,1)  = gretna_auc(Net_para.bin.ass(:,isub),Delta);
    Net_para.bin.abeta(isub,1) = gretna_auc(Net_para.bin.beta(:,isub),Delta);
    Net_para.bin.aq(isub,1)    = gretna_auc(Net_para.bin.q(:,isub),Delta);
    
    Net_para.bin.acpratio(isub,1)   = gretna_auc(Net_para.bin.cpratio(:,isub),Delta);
    Net_para.bin.alpratio(isub,1)   = gretna_auc(Net_para.bin.lpratio(:,isub),Delta);
    Net_para.bin.aloceratio(isub,1) = gretna_auc(Net_para.bin.loceratio(:,isub),Delta);
    Net_para.bin.ageratio(isub,1)   = gretna_auc(Net_para.bin.geratio(:,isub),Delta);
    Net_para.bin.adegratio(isub,1)  = gretna_auc(Net_para.bin.degratio(:,isub),Delta);
    Net_para.bin.abwratio(isub,1)   = gretna_auc(Net_para.bin.bwratio(:,isub),Delta);
    Net_para.bin.abgratio(isub,1)   = gretna_auc(Net_para.bin.bgratio(:,isub),Delta);
    Net_para.bin.acloratio(isub,1)  = gretna_auc(Net_para.bin.cloratio(:,isub),Delta);
    Net_para.bin.acomratio(isub,1)  = gretna_auc(Net_para.bin.comratio(:,isub),Delta);
    Net_para.bin.aeigratio(isub,1)  = gretna_auc(Net_para.bin.eigratio(:,isub),Delta);
    Net_para.bin.ahindratio(isub,1) = gretna_auc(Net_para.bin.hindratio(:,isub),Delta);
    Net_para.bin.ainfratio(isub,1)  = gretna_auc(Net_para.bin.infratio(:,isub),Delta);
    Net_para.bin.akatzratio(isub,1) = gretna_auc(Net_para.bin.katzratio(:,isub),Delta);
    Net_para.bin.alapratio(isub,1)  = gretna_auc(Net_para.bin.lapratio(:,isub),Delta);
    Net_para.bin.alevratio(isub,1)  = gretna_auc(Net_para.bin.levratio(:,isub),Delta);
    Net_para.bin.apageratio(isub,1) = gretna_auc(Net_para.bin.pageratio(:,isub),Delta);
    Net_para.bin.aassratio(isub,1)  = gretna_auc(Net_para.bin.assratio(:,isub),Delta);
    Net_para.bin.abetaratio(isub,1) = gretna_auc(Net_para.bin.betaratio(:,isub),Delta);
    Net_para.bin.aqratio(isub,1)    = gretna_auc(Net_para.bin.qratio(:,isub),Delta);
    
    
    Net_para.wei.acp(isub,1)   = gretna_auc(Net_para.wei.cp(:,isub),Delta);
    Net_para.wei.alp(isub,1)   = gretna_auc(Net_para.wei.lp(:,isub),Delta);
    Net_para.wei.aloce(isub,1) = gretna_auc(Net_para.wei.loce(:,isub),Delta);
    Net_para.wei.age(isub,1)   = gretna_auc(Net_para.wei.ge(:,isub),Delta);
    Net_para.wei.adeg(isub,1)  = gretna_auc(Net_para.wei.deg(:,isub),Delta);
    Net_para.wei.abw(isub,1)   = gretna_auc(Net_para.wei.bw(:,isub),Delta);
    Net_para.wei.abg(isub,1)   = gretna_auc(Net_para.wei.bg(:,isub),Delta);
    Net_para.wei.aclo(isub,1)  = gretna_auc(Net_para.wei.clo(:,isub),Delta);
    Net_para.wei.acom(isub,1)  = gretna_auc(Net_para.wei.com(:,isub),Delta);
    Net_para.wei.aeig(isub,1)  = gretna_auc(Net_para.wei.eig(:,isub),Delta);
    Net_para.wei.ahind(isub,1) = gretna_auc(Net_para.wei.hind(:,isub),Delta);
    Net_para.wei.ainf(isub,1)  = gretna_auc(Net_para.wei.inf(:,isub),Delta);
    Net_para.wei.akatz(isub,1) = gretna_auc(Net_para.wei.katz(:,isub),Delta);
    Net_para.wei.alap(isub,1)  = gretna_auc(Net_para.wei.lap(:,isub),Delta);
    Net_para.wei.alev(isub,1)  = gretna_auc(Net_para.wei.lev(:,isub),Delta);
    Net_para.wei.apage(isub,1) = gretna_auc(Net_para.wei.page(:,isub),Delta);
    Net_para.wei.aass(isub,1)  = gretna_auc(Net_para.wei.ass(:,isub),Delta);
    Net_para.wei.abeta(isub,1) = gretna_auc(Net_para.wei.beta(:,isub),Delta);
    Net_para.wei.aq(isub,1)    = gretna_auc(Net_para.wei.q(:,isub),Delta);
    
    Net_para.wei.acpratio(isub,1)   = gretna_auc(Net_para.wei.cpratio(:,isub),Delta);
    Net_para.wei.alpratio(isub,1)   = gretna_auc(Net_para.wei.lpratio(:,isub),Delta);
    Net_para.wei.aloceratio(isub,1) = gretna_auc(Net_para.wei.loceratio(:,isub),Delta);
    Net_para.wei.ageratio(isub,1)   = gretna_auc(Net_para.wei.geratio(:,isub),Delta);
    Net_para.wei.adegratio(isub,1)  = gretna_auc(Net_para.wei.degratio(:,isub),Delta);
    Net_para.wei.abwratio(isub,1)   = gretna_auc(Net_para.wei.bwratio(:,isub),Delta);
    Net_para.wei.abgratio(isub,1)   = gretna_auc(Net_para.wei.bgratio(:,isub),Delta);
    Net_para.wei.acloratio(isub,1)  = gretna_auc(Net_para.wei.cloratio(:,isub),Delta);
    Net_para.wei.acomratio(isub,1)  = gretna_auc(Net_para.wei.comratio(:,isub),Delta);
    Net_para.wei.aeigratio(isub,1)  = gretna_auc(Net_para.wei.eigratio(:,isub),Delta);
    Net_para.wei.ahindratio(isub,1) = gretna_auc(Net_para.wei.hindratio(:,isub),Delta);
    Net_para.wei.ainfratio(isub,1)  = gretna_auc(Net_para.wei.infratio(:,isub),Delta);
    Net_para.wei.akatzratio(isub,1) = gretna_auc(Net_para.wei.katzratio(:,isub),Delta);
    Net_para.wei.alapratio(isub,1)  = gretna_auc(Net_para.wei.lapratio(:,isub),Delta);
    Net_para.wei.alevratio(isub,1)  = gretna_auc(Net_para.wei.levratio(:,isub),Delta);
    Net_para.wei.apageratio(isub,1) = gretna_auc(Net_para.wei.pageratio(:,isub),Delta);
    Net_para.wei.aassratio(isub,1)  = gretna_auc(Net_para.wei.assratio(:,isub),Delta);
    Net_para.wei.abetaratio(isub,1) = gretna_auc(Net_para.wei.betaratio(:,isub),Delta);
    Net_para.wei.aqratio(isub,1)    = gretna_auc(Net_para.wei.qratio(:,isub),Delta);
    
    Node_para.bin.acp(isub,:)   = gretna_auc(squeeze(Node_para.bin.cp(:,isub,:)),Delta);
    Node_para.bin.alp(isub,:)   = gretna_auc(squeeze(Node_para.bin.lp(:,isub,:)),Delta);
    Node_para.bin.aloce(isub,:) = gretna_auc(squeeze(Node_para.bin.loce(:,isub,:)),Delta);
    Node_para.bin.age(isub,:)   = gretna_auc(squeeze(Node_para.bin.ge(:,isub,:)),Delta);
    Node_para.bin.adeg(isub,:)  = gretna_auc(squeeze(Node_para.bin.deg(:,isub,:)),Delta);
    Node_para.bin.abw(isub,:)   = gretna_auc(squeeze(Node_para.bin.bw(:,isub,:)),Delta);
    Node_para.bin.abg(isub,:)   = gretna_auc(squeeze(Node_para.bin.bg(:,isub,:)),Delta);
    Node_para.bin.aclo(isub,:)  = gretna_auc(squeeze(Node_para.bin.clo(:,isub,:)),Delta);
    Node_para.bin.acom(isub,:)  = gretna_auc(squeeze(Node_para.bin.com(:,isub,:)),Delta);
    Node_para.bin.aeig(isub,:)  = gretna_auc(squeeze(Node_para.bin.eig(:,isub,:)),Delta);
    Node_para.bin.ahind(isub,:) = gretna_auc(squeeze(Node_para.bin.hind(:,isub,:)),Delta);
    Node_para.bin.ainf(isub,:)  = gretna_auc(squeeze(Node_para.bin.inf(:,isub,:)),Delta);
    Node_para.bin.akatz(isub,:) = gretna_auc(squeeze(Node_para.bin.katz(:,isub,:)),Delta);
    Node_para.bin.alap(isub,:)  = gretna_auc(squeeze(Node_para.bin.lap(:,isub,:)),Delta);
    Node_para.bin.alev(isub,:)  = gretna_auc(squeeze(Node_para.bin.lev(:,isub,:)),Delta);
    Node_para.bin.apage(isub,:) = gretna_auc(squeeze(Node_para.bin.page(:,isub,:)),Delta);
    
    Node_para.bin.acpratio(isub,:)   = gretna_auc(squeeze(Node_para.bin.cpratio(:,isub,:)),Delta);
    Node_para.bin.alpratio(isub,:)   = gretna_auc(squeeze(Node_para.bin.lpratio(:,isub,:)),Delta);
    Node_para.bin.aloceratio(isub,:) = gretna_auc(squeeze(Node_para.bin.loceratio(:,isub,:)),Delta);
    Node_para.bin.ageratio(isub,:)   = gretna_auc(squeeze(Node_para.bin.geratio(:,isub,:)),Delta);
    Node_para.bin.adegratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.degratio(:,isub,:)),Delta);
    Node_para.bin.abwratio(isub,:)   = gretna_auc(squeeze(Node_para.bin.bwratio(:,isub,:)),Delta);
    Node_para.bin.abgratio(isub,:)   = gretna_auc(squeeze(Node_para.bin.bgratio(:,isub,:)),Delta);
    Node_para.bin.acloratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.cloratio(:,isub,:)),Delta);
    Node_para.bin.acomratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.comratio(:,isub,:)),Delta);
    Node_para.bin.aeigratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.eigratio(:,isub,:)),Delta);
    Node_para.bin.ahindratio(isub,:) = gretna_auc(squeeze(Node_para.bin.hindratio(:,isub,:)),Delta);
    Node_para.bin.ainfratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.infratio(:,isub,:)),Delta);
    Node_para.bin.akatzratio(isub,:) = gretna_auc(squeeze(Node_para.bin.katzratio(:,isub,:)),Delta);
    Node_para.bin.alapratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.lapratio(:,isub,:)),Delta);
    Node_para.bin.alevratio(isub,:)  = gretna_auc(squeeze(Node_para.bin.levratio(:,isub,:)),Delta);
    Node_para.bin.apageratio(isub,:) = gretna_auc(squeeze(Node_para.bin.pageratio(:,isub,:)),Delta);
    
    
    Node_para.wei.acp(isub,:)   = gretna_auc(squeeze(Node_para.wei.cp(:,isub,:)),Delta);
    Node_para.wei.alp(isub,:)   = gretna_auc(squeeze(Node_para.wei.lp(:,isub,:)),Delta);
    Node_para.wei.aloce(isub,:) = gretna_auc(squeeze(Node_para.wei.loce(:,isub,:)),Delta);
    Node_para.wei.age(isub,:)   = gretna_auc(squeeze(Node_para.wei.ge(:,isub,:)),Delta);
    Node_para.wei.adeg(isub,:)  = gretna_auc(squeeze(Node_para.wei.deg(:,isub,:)),Delta);
    Node_para.wei.abw(isub,:)   = gretna_auc(squeeze(Node_para.wei.bw(:,isub,:)),Delta);
    Node_para.wei.abg(isub,:)   = gretna_auc(squeeze(Node_para.wei.bg(:,isub,:)),Delta);
    Node_para.wei.aclo(isub,:)  = gretna_auc(squeeze(Node_para.wei.clo(:,isub,:)),Delta);
    Node_para.wei.acom(isub,:)  = gretna_auc(squeeze(Node_para.wei.com(:,isub,:)),Delta);
    Node_para.wei.aeig(isub,:)  = gretna_auc(squeeze(Node_para.wei.eig(:,isub,:)),Delta);
    Node_para.wei.ahind(isub,:) = gretna_auc(squeeze(Node_para.wei.hind(:,isub,:)),Delta);
    Node_para.wei.ainf(isub,:)  = gretna_auc(squeeze(Node_para.wei.inf(:,isub,:)),Delta);
    Node_para.wei.akatz(isub,:) = gretna_auc(squeeze(Node_para.wei.katz(:,isub,:)),Delta);
    Node_para.wei.alap(isub,:)  = gretna_auc(squeeze(Node_para.wei.lap(:,isub,:)),Delta);
    Node_para.wei.alev(isub,:)  = gretna_auc(squeeze(Node_para.wei.lev(:,isub,:)),Delta);
    Node_para.wei.apage(isub,:) = gretna_auc(squeeze(Node_para.wei.page(:,isub,:)),Delta);
    
    Node_para.wei.acpratio(isub,:)   = gretna_auc(squeeze(Node_para.wei.cpratio(:,isub,:)),Delta);
    Node_para.wei.alpratio(isub,:)   = gretna_auc(squeeze(Node_para.wei.lpratio(:,isub,:)),Delta);
    Node_para.wei.aloceratio(isub,:) = gretna_auc(squeeze(Node_para.wei.loceratio(:,isub,:)),Delta);
    Node_para.wei.ageratio(isub,:)   = gretna_auc(squeeze(Node_para.wei.geratio(:,isub,:)),Delta);
    Node_para.wei.adegratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.degratio(:,isub,:)),Delta);
    Node_para.wei.abwratio(isub,:)   = gretna_auc(squeeze(Node_para.wei.bwratio(:,isub,:)),Delta);
    Node_para.wei.abgratio(isub,:)   = gretna_auc(squeeze(Node_para.wei.bgratio(:,isub,:)),Delta);
    Node_para.wei.acloratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.cloratio(:,isub,:)),Delta);
    Node_para.wei.acomratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.comratio(:,isub,:)),Delta);
    Node_para.wei.aeigratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.eigratio(:,isub,:)),Delta);
    Node_para.wei.ahindratio(isub,:) = gretna_auc(squeeze(Node_para.wei.hindratio(:,isub,:)),Delta);
    Node_para.wei.ainfratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.infratio(:,isub,:)),Delta);
    Node_para.wei.akatzratio(isub,:) = gretna_auc(squeeze(Node_para.wei.katzratio(:,isub,:)),Delta);
    Node_para.wei.alapratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.lapratio(:,isub,:)),Delta);
    Node_para.wei.alevratio(isub,:)  = gretna_auc(squeeze(Node_para.wei.levratio(:,isub,:)),Delta);
    Node_para.wei.apageratio(isub,:) = gretna_auc(squeeze(Node_para.wei.pageratio(:,isub,:)),Delta);
end

return