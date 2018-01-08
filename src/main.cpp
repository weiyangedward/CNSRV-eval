/*
 1. This code is used to test the noise cluster in a clustering results.
 Specifically, it aims to reason why the size of noise cluster is very
 small in the clustering result (noise size < 10).
 
 2. To justify this result, this code takes a node out of a cluster at
 the time, and see the change of final score, if the score decreases,
 then this node should not be put into the currently cluster and should
 be in the nosie cluster instead. If the clustering algorithm works
 correctly, this case should be very rare.
 
 3. Note that this code is not suitable for the new InOutRatio clustering
 algorithm, where the objective function is changed to:
 
 score =  sum_over_i(in_cluster_density(i) / ( in_cluster_density(i) + out_cluster_density(i))) * n + lamda * noise_term
 
 where lamda = sum_over_i(in_cluster_density(i) / ( in_cluster_density(i) + out_cluster_density(i))) / N
 
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <list>
#include <unordered_set>
#include <stdio.h>
using namespace std;

void inOutCost(vector<vector<unordered_set<int>>> & spe_clus, vector<unordered_map<int, vector<int>>> & adj_list, vector<int> & total_gene_num, char* output_file)
{
    double ave_in_density, ave_out_density, ave_in_out_density_ratio;
    int count = 0;
    FILE * pFile2;
    pFile2 = fopen(output_file,"w");
    
    for (int i = 0; i<3;i++) // 3 spe
    {
        cout << "spe " << i << endl;
        for (int j=0;j<=10;j++) // 11 clusters
        {
            if (spe_clus[i][j].size() > 0) // if number of nodes > 0 in spe i cluster j
            {
                /*========== compute in/out density of cluster ============*/
                int in_edge = 0, out_edge = 0, num_node = (int)spe_clus[i][j].size();
                for (auto cur_gene:spe_clus[i][j])
                {
                    for (auto gene2:adj_list[i].at(cur_gene)) // for each gene connects to cur_gene
                    {
                        if (spe_clus[i][j].find(gene2) != spe_clus[i][j].end()) // if gene2 in cluster [i][j]
                            in_edge++;
                        else
                            out_edge++;
                    }
                }
                in_edge /= 2; // in_edge reduce to half because it was double counted in cluster
                int possible_in_edge = num_node * (num_node-1) / 2;
                int possible_out_edge = num_node * (total_gene_num[i] - num_node);
                double in_density = (double)in_edge / (num_node * (num_node-1) / 2);
                
                double out_density = (double)out_edge / (num_node * (total_gene_num[i] - num_node));
                double in_out_density_ratio = in_density / out_density;
                double in_vs_inPlusOut_density_ratio =  in_density / (in_density + out_density);
                double cost = ( in_density / (in_density + out_density) ) * num_node;
                
                fprintf(pFile2, "spe \t %d \t cluster \t %d \t #nodes \t %d \t in-edge \t %d \t possible_in_edge \t %d \t in-density \t %g \t out-edge \t %d \t possible_out_edge \t %d \t out-density \t %g \t in/out-density ratio \t %g \t in/(in+out) density ratio \t %g \t cost \t %g\n", i, j, num_node, in_edge,possible_in_edge, in_density, out_edge,possible_out_edge, out_density,in_out_density_ratio, in_vs_inPlusOut_density_ratio, cost);
                if (j != 0)
                {
                    count ++;
                    ave_in_density += in_density;
                    ave_out_density += out_density;
                    ave_in_out_density_ratio += in_out_density_ratio;
                }
            }
        }
    }
    ave_in_density /= count;
    ave_out_density /= count;
    ave_in_out_density_ratio /= count;
    fprintf(pFile2, "average_excluding_noise_cluster\t ave_in_density \t %g \t ave_out_density \t %g \t ave_in_out_density_ratio \t %g\n", ave_in_density, ave_out_density, ave_in_out_density_ratio);
    
    fclose(pFile2);
}


int main(int argc, char * argv[]) {
    if (argc != 4){
        printf("usage:\n\t%s <clustering output> <network> <outputFile>\n", argv[0]);
        exit(1);
    }
    
    /*=========== read clustering output file ===========*/
    ifstream CLUSTEROUT(argv[1]);
    if (!CLUSTEROUT){
        cerr << "clustering output could not be opened for reading!" << endl;
        exit(1);
    }
    
    string title; // "Best clustering..."
    getline(CLUSTEROUT, title);
    
    vector<vector<unordered_set<int>>> spe_clus; // arr to store genes per cluster per species
    vector<int> total_gene_num(3,0); // arr to store the total number of genes per species
    
    // initialize spe_clus
    for (int i = 0; i<3;i++)
    {
        vector<unordered_set<int>> spe;
        for (int j = 0; j<=10; j++)
        {
            unordered_set<int> clus;
            spe.push_back(clus);
        }
        spe_clus.push_back(spe);
    }
    
    // read file CLUSTEROUT
    while (CLUSTEROUT)
    {
        string line;
        getline(CLUSTEROUT, line);
        istringstream split_line(line);
        string w1, w2, w3;
        int spe_id, gene_id, cluster_id;
        split_line >> w1 >> spe_id >> w2 >> gene_id >> w3 >> cluster_id;
        
        total_gene_num[spe_id]++;
        
        spe_clus[spe_id][cluster_id].emplace(gene_id);
    }
    CLUSTEROUT.close();
    printf("total number of genes: hb:%d, mm:%d, sb:%d\n", total_gene_num[0], total_gene_num[1], total_gene_num[2]);
    
    
    /*=========== read network file ==============*/
    ifstream NETWORK(argv[2]);
    if (!NETWORK){
        cerr << "network file could not be opened for reading!" << endl;
        exit(1);
    }
    
    getline(NETWORK, title); // "Real network"
    
    unordered_map<int, vector<int>> spe0_adj_list;
    unordered_map<int, vector<int>> spe1_adj_list;
    unordered_map<int, vector<int>> spe2_adj_list;
    vector<unordered_map<int, vector<int>>> adj_list; // arr to store hash map for adjacency list per gene per spe
    adj_list.push_back(spe0_adj_list); adj_list.push_back(spe1_adj_list); adj_list.push_back(spe2_adj_list);
    
    // read file NETWORK
    while (NETWORK)
    {
        string line;
        getline(NETWORK, line);
        istringstream split_line(line);
        int spe, gene_id1;
        int gene_id2, edge;
        split_line >> spe >> gene_id1 >> gene_id2 >> edge;
        
        if (adj_list[spe].find(gene_id1) == adj_list[spe].end()) // if key does not exist in hash table
        {
            vector<int> tmp;
            tmp.push_back(gene_id2);
            adj_list[spe].emplace(gene_id1, tmp);
        }
        else // if key exists
        {
            adj_list[spe].at(gene_id1).push_back(gene_id2);
        }
    }
    NETWORK.close();
    
    /*========= compute in-out cost per cluster =========*/
    inOutCost(spe_clus, adj_list, total_gene_num, argv[3]);
    
    return 0;
}
