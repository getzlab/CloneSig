import pandas as pd
import numpy as np
import os
import sys
import scipy
import scipy.stats as stat
import scipy.sparse as sp
import numpy_groupies as npg
import multiprocessing
sys.path.append("CurveBall/cbsig")
import curveball

class CloneSigAnalyzer:
    def __init__(self, input_maf, gene_list_txt):
        """
        Initialize the CloneigAnalyzer with input data
        
        Args:
            input_maf (str): Path to the MAF file
            gene_list_txt (str): Path to the gene list file
        """
        self.input_maf = input_maf
        self.gene_list = pd.read_csv(gene_list_txt, sep='\t')['genes'].tolist()
        self.maf = pd.read_csv(input_maf, sep='\t')
        self.maf = self.maf.sort_values(by=['Tumor_Sample_Barcode']).reset_index(drop=True)
        self.subclones_sig_genes = pd.read_csv('sig_genes_subclones.txt', sep='\t')
        
        # Setup basic matrices and indices
        self._setup_matrices()
        
        # Default parameters
        self.n_perm = 100000
        self.permutation_dir = "permutation_matrices"
        self.batch_size = 5

    def _setup_matrices(self):
        """Setup the basic matrices and indices needed for analysis"""
        clones, c_ui, c_idx = np.unique(self.maf["Tumor_Sample_Barcode"], return_inverse=True, return_index=True)
        _, p_idx = np.unique(self.maf["Patient_ID"], return_inverse=True)
        self.c2pat = p_idx[c_ui]
        genes, g_ui, g_idx = np.unique(self.maf["Hugo_Symbol"], return_inverse=True, return_index=True)
        self.g2gidx = dict(zip(genes, g_idx[g_ui]))
        
        # Create initial CG matrix
        self.CG = np.array((sp.coo_matrix((np.ones_like(c_idx), (c_idx, g_idx))).todense() > 0))
        self.hp = curveball.find_presences(self.CG)

    def calculate_distributions(self, gene):
        """
        Calculate probability distributions for mutations in a gene
        
        Args:
            gene (str): Gene of interest
            
        Returns:
            DataFrame with probability distributions
        """
        input_maf = self.maf.copy()
        all_patients = list(input_maf['Patient_ID'].unique())
        
        # Total number of mutations in the target gene is N0
        gene_all_muts = input_maf[input_maf['Hugo_Symbol']==gene]
        N0 = len(gene_all_muts.drop_duplicates(['Tumor_Sample_Barcode']))
        total_muts = len(input_maf)
        all_probs = []
        
        for patient in all_patients: 
            RA_maf = input_maf[input_maf['Patient_ID']==patient]
            RA_maf_clones = list(RA_maf['Tumor_Sample_Barcode'].unique())
            P_pat = np.zeros((len(RA_maf_clones), 2))
            
            for clone_idx, clone in enumerate(RA_maf_clones):
                clone_df = RA_maf[RA_maf['Tumor_Sample_Barcode'] == clone]
                Mij = len(clone_df)
                term_to_calc = (1 - (Mij/total_muts))**N0
                
                P_pat[clone_idx, 0] = term_to_calc
                P_pat[clone_idx, 1] = 1-term_to_calc

            P_prob_df = pd.DataFrame(P_pat, columns=[0, 1])
            P_prob_df.insert(0, 'clone', RA_maf_clones)
            all_probs.append(P_prob_df)
        
        all_probs_concat = pd.concat(all_probs)
        pats_label = [clone.split('_')[0] for clone in all_probs_concat['clone']]
        
        all_probs_df_labeled = all_probs_concat.copy()
        all_probs_df_labeled['patient'] = pats_label
        
        clone_order = list(np.sort(input_maf["Tumor_Sample_Barcode"]))
        order_map = {val:i for i,val in enumerate(clone_order)}
        all_probs_df_labeled['clones_sort'] = all_probs_df_labeled['clone'].map(order_map)
        all_probs_df_labeled = all_probs_df_labeled.sort_values(by='clones_sort')
        all_probs_df_labeled = all_probs_df_labeled.drop('clones_sort', axis=1)
        all_probs_df_labeled.reset_index(inplace=True, drop=True)
        
        return all_probs_df_labeled

    def calculate_scores_all_identity_vectors(self, mutated_clone_per_pat, gene_pm):
        """
        Calculate mutation scores per patient
        
        Args:
            gene_pm: DataFrame with probability distribution of mutations per clone
            mutated_clone_per_pat: Dict mapping patient index to list of mutated clone indices
            
        Returns:
            Array of scores
        """
        unique_pats = gene_pm['patient'].unique()
        pat_dict = {pat: idx for idx, pat in enumerate(unique_pats)}
        score_vector = np.zeros(len(unique_pats))
        
        log_prob_0 = np.log(gene_pm[0])
        log_prob_1 = np.log(gene_pm[1])
        
        for pat in unique_pats:
            pat_mask = gene_pm['patient'] == pat
            pat_probs = gene_pm[pat_mask]
            pat_num = pat_dict[pat]
            count_muts = len(mutated_clone_per_pat[pat_num])
            n_clones = len(pat_probs)

            if count_muts <= 1:
                log_prob_1_masked = log_prob_1.values[pat_mask]
                log_prob_0_masked = log_prob_0.values[pat_mask]
                total_score = np.sum(log_prob_0_masked)
                
                identity_matrix = np.eye(n_clones, dtype=int)
                single_mut_scores = np.sum(log_prob_1_masked * identity_matrix, axis=1) + \
                                  np.sum(log_prob_0_masked * (1 - identity_matrix), axis=1)
                
                all_null_scores = np.concatenate([single_mut_scores, [total_score]])
                score_vector[pat_num] = -scipy.special.logsumexp(all_null_scores)
            else:
                mut_clones_idx = set(mutated_clone_per_pat[pat_num])
                pat_indices = pat_probs.index
                is_mutated = np.array([idx in mut_clones_idx for idx in pat_indices])
                score_vector[pat_num] = -np.sum(log_prob_1.values[pat_mask][is_mutated]) - \
                                      np.sum(log_prob_0.values[pat_mask][~is_mutated])

        return score_vector

    def compute_mut_clone_indices(self, pat_idx, mat):
        """
        Compute mutation clone indices
        
        Args:
            pat_idx: Patient indices
            mat: Mutation matrix
            
        Returns:
            Dictionary of mutation clone indices per patient
        """
        mut_clone_indices = {pat: [] for pat in np.unique(pat_idx)}
        
        for clone_idx, pat in enumerate(pat_idx):
            mutated_clones = np.where(mat[clone_idx, :] > 0)[0]
            mut_clone_indices[pat].extend([clone_idx for gene in mutated_clones])
            
        return mut_clone_indices

    def process_single_gene(self, gene):
        """
        Process a single gene and return its results
        
        Args:
            gene: Gene to process
            
        Returns:
            Dictionary with gene results
        """
        print(f"Processing gene: {gene}")
        
        # Calculate observed distribution
        gene_Pm = self.calculate_distributions(gene)
        
        # Get gene indices
        subset_idx = (self.subclones_sig_genes['gene'] == gene)
        subset_idx_mat = [self.g2gidx[x] for x in self.subclones_sig_genes.loc[subset_idx, "gene"]]
        
        # Calculate observed score
        mut_clone_indices_obs = self.compute_mut_clone_indices(self.c2pat, self.CG[:, subset_idx_mat])
        obs_score_array = self.calculate_scores_all_identity_vectors(mut_clone_indices_obs, gene_Pm)
        obs_score = np.sum(obs_score_array)
        
        # Process permutations
        test = 0
        test_exact = 0
        n_perm = 0
        
        for chunk_file in sorted(os.listdir(self.permutation_dir)):
            if not chunk_file.endswith('.npz'):
                continue
                
            chunk_path = os.path.join(self.permutation_dir, chunk_file)
            chunk_matrices = np.load(chunk_path)['matrices']
            
            for perm_matrix in chunk_matrices:
                itr_test = self.compute_mut_clone_indices(self.c2pat, perm_matrix[:, subset_idx_mat])
                score_per_itr_per_pat = self.calculate_scores_all_identity_vectors(itr_test, gene_Pm)
                score_per_itr = np.sum(score_per_itr_per_pat)
                
                test += (score_per_itr >= obs_score)
                test_exact += (np.abs(score_per_itr) == np.abs(obs_score))
                n_perm += 1
        
        p_values = test / n_perm
        p_val_exact = test_exact / n_perm
        
        return {
            'gene': gene,
            'Cb_p_vals': float(p_values),
            'Cb_exact_p_vals': float(p_val_exact)
        }

    def generate_chunk_permutations(self, chunk_id, chunk_size):
        """
        Generate a chunk of permutations
        
        Args:
            chunk_id: ID of the chunk
            chunk_size: Size of the chunk
            
        Returns:
            chunk_id
        """
        print(f"Processing chunk {chunk_id}")
        chunk_matrices = []
        
        # Do a longer burn-in at the start of chunk
        CG_current = curveball.curve_ball(self.CG, self.hp, 10000)
        hp_current = curveball.find_presences(CG_current)
        
        # Generate permutations for this chunk
        for i in range(chunk_size):
            CG_current = curveball.curve_ball(CG_current, hp_current, 10000)
            chunk_matrices.append(CG_current)
            
            if (i + 1) % 100 == 0:
                print(f"Chunk {chunk_id}: Generated {i + 1}/{chunk_size} permutations")
        
        # Save chunk to disk
        chunk_file = os.path.join(self.permutation_dir, f"permutations_chunk_{chunk_id}.npz")
        np.savez_compressed(chunk_file, matrices=chunk_matrices)
        del chunk_matrices
        return chunk_id

    def generate_permutations(self):
        """Generate all permutations in parallel"""
        os.makedirs(self.permutation_dir, exist_ok=True)
        
        # Initial burn-in
        CG_burnin = curveball.curve_ball(self.CG, self.hp, 10000)
        hp_burnin = curveball.find_presences(CG_burnin)
        
        chunk_size = 1000
        n_chunks = self.n_perm // chunk_size
        
        n_processes = min(n_chunks, max(1, multiprocessing.cpu_count() // 2))
        print(f"Generating permutations using {n_processes} processes...")
        
        with multiprocessing.Pool(processes=n_processes) as pool:
            results = pool.starmap(
                self.generate_chunk_permutations,
                [(chunk_id, chunk_size) for chunk_id in range(n_chunks)]
            )
        
        print(f"Completed generating {n_chunks} chunks of permutations")

    def process_genes_parallel(self):
        """
        Process all genes in parallel batches
        
        Returns:
            List of results for all genes
        """
        results = []
        for i in range(0, len(self.gene_list), self.batch_size):
            batch = self.gene_list[i:i + self.batch_size]
            with multiprocessing.Pool(processes=min(len(batch), multiprocessing.cpu_count())) as pool:
                batch_results = pool.map(self.process_single_gene, batch)
            results.extend(batch_results)
        return results

    def run_analysis(self, step='all'):
        """
        Run the complete analysis
        
        Args:
            step: Which step to run ('permutations', 'genes', or 'all')
        """
        if step == 'all' or step == 'permutations':
            self.generate_permutations()
            
        if step == 'all' or step == 'genes':
            gene_results = self.process_genes_parallel()
            Cb_all_genes_df = pd.DataFrame(gene_results)
            Cb_all_genes_df.to_csv('Cb_onesided_subclones_040325_100k.txt', sep='\t', index=False) 