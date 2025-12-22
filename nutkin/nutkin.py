import os
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import itertools
from sklearn.decomposition import PCA
import seaborn as sns
from config import *
import anndata as ad

class Nutkin:
    def __init__(self, adata, output_dir, num_pc = DEFAULT_NUM_PC, group_col="group", metric_type="sum",verbal = True,detail=False):
        self.verbal = verbal
        self.detail = detail
        self.group_col = group_col
        self.output_dir = output_dir
        self.num_pc = num_pc
        self.metric_type = metric_type.lower()
        if self.metric_type not in ["sum", "product", "both"]:
            raise ValueError("metric_type must be 'sum', 'product', or 'both'")
            
            
        os.makedirs(self.output_dir, exist_ok=True)
        if group_col not in adata.obs.columns:
            raise ValueError(f"'{group_col}' not found in adata.obs.")
        self.adata = adata
        self.data = self.adata_to_df_with_obs(adata, group_col)
        self.dimReduceData = pd.DataFrame()
        
    def adata_to_df_with_obs(self,adata, group_col):
    
        if isinstance(group_col, str):
            group_col = [group_col]
    
        # Convert expression matrix to DataFrame
        if hasattr(adata.X, "toarray"):  # sparse matrix case
            expr_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
        else:
            expr_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    
        # Select obs columns
        obs_df = adata.obs[group_col]
    
        # Concatenate expression and obs info
        merged_df = pd.concat([obs_df.reset_index(drop=True), expr_df.reset_index(drop=True)], axis=1)
    
        return merged_df
    
    def _optional_print(self,msg,*msgs):
        if self.verbal:
            additional_msg = ""
            for m in msgs:
                additional_msg += str(m) 
            combined_msg = msg + additional_msg
            print(combined_msg+'\n')
        return

    def _dim_reduction(self, method="pca"):
        if not self.dimReduceData.empty:
            return self.dimReduceData
        
        self._optional_print("Running dimensional reduction (PCA)...")
        n_components = min(self.num_pc, self.data.shape[1] - 1)
        self.transformer = PCA(n_components)
        
        dimReduceData = pd.DataFrame(self.transformer.fit_transform(self.data.iloc[:, 1:]))
        self.dimReduceData = pd.concat([self.data.iloc[:,0],dimReduceData], axis=1)
        
        return self.dimReduceData

   
    def _change_of_axes(self, group_data):
        """represent group of cell into group specific features using pca"""
        transformed_group_data = pd.DataFrame(self.transformer.fit_transform(group_data.iloc[:, 1:]))
        return pd.concat([group_data.iloc[:,0].reset_index(drop=True),transformed_group_data.reset_index(drop=True)], axis=1).reset_index(drop=True)
    
    
    def _summarize_cov_matrix(self, group_matrix):
        """Calculate variability metrics based on metric_type"""
        var_matrix = np.cov(group_matrix.T)
        
        if isinstance(var_matrix, float):
            diag = np.array([var_matrix])
        else:
            diag = np.diagonal(var_matrix).copy()
        
        eps = 1e-50
        diag[diag < eps] = eps
        
        metric_sum = np.log(np.sum(np.sqrt(diag)))
        metric_product = np.sum(np.log(np.sqrt(diag)))
        
        if self.metric_type == "sum":
            return metric_sum
        elif self.metric_type == "product":
            return metric_product
        else:  # both
            return metric_sum, metric_product
    
    
    def measure(self):
        self._optional_print("Quantifying variability")
        self._dim_reduction()
        output_df = pd.DataFrame()
        groups = self.data.iloc[:,0].unique()
  
        self.MIN_GROUP_SIZE = 3
    
        for group in groups:
            group_cells = self.data[self.data.iloc[:, 0] == group]
            num_cell = len(group_cells)

        
            if num_cell < self.MIN_GROUP_SIZE:
                self._optional_print(f"Skipping group '{group}' (too few cells: {num_cell})")
                continue
            metric, num_cell = self._measure_one_group(group)
            
            if self.metric_type == "both":
                metric_sum, metric_product = metric
                sort_col = "sum statistics"
                metric_df = pd.DataFrame([{
                    "group": group, 
                    "number of cells": num_cell, 
                    "sum statistics": round(metric_sum, DECIMAL),
                    "product statistics": round(metric_product, DECIMAL)
                }])
            else:
                sort_col = f"{self.metric_type} statistics"
                metric_df = pd.DataFrame([{
                    "group": group,
                    "number of cells": num_cell,
                    f"{self.metric_type} statistics": round(metric, DECIMAL)
                }])
            output_df = pd.concat([output_df, metric_df], ignore_index=True)
        output_df = output_df.sort_values(by=sort_col, ascending=True).reset_index(drop=True)
 
        self._optional_print("Saving results")
        output_df.to_csv(os.path.join(self.output_dir, "var_measurement.csv"),index=False)
        return output_df
    
    def _measure_one_group(self,group):
        if self.dimReduceData.empty == True:
            self._dim_reduction()
        group_data = self.dimReduceData.loc[self.dimReduceData[self.dimReduceData.columns[0]] == group]
        group_matrix = self._change_of_axes(group_data)
        group_size=group_data.shape[0]
        metric = self._summarize_cov_matrix(group_matrix.iloc[:, 1:])
        return metric,group_size

    def _run_bootstrap(self, group1,group2,num_resample,seed = 42):
        self._optional_print("Running Bootstrapping")
        random.seed(seed)
        self._dim_reduction()
        label_col_name = self.dimReduceData.columns.tolist()[0]
        group1_data = self.dimReduceData.loc[self.dimReduceData[label_col_name]==group1]
        group2_data = self.dimReduceData.loc[self.dimReduceData[label_col_name]==group2]
        
        num_resample_group2=num_resample
        num_resample_group1=num_resample
        
        
        g1_sum_stats, g2_sum_stats = [], []
        g1_product_stats, g2_product_stats = [], []
        optimal_bootstrap_sample_size= min(BOOTSTRAP_SAMPLE_SIZE_UPPER_BOUND, int(MAXIMUM_BOOTSTRAP_PROPORTION*min(len(group1_data), len(group2_data))))
        for rep in range(num_resample):
            resample1 = group1_data.sample(optimal_bootstrap_sample_size, replace=True, random_state=seed + rep)
            resample2 = group2_data.sample(optimal_bootstrap_sample_size, replace=True, random_state=seed + rep)

            changed1 = self._change_of_axes(resample1)
            changed2 = self._change_of_axes(resample2)

            metric1 = self._summarize_cov_matrix(changed1.iloc[:, 1:])
            metric2 = self._summarize_cov_matrix(changed2.iloc[:, 1:])

            if self.metric_type == "both":
                g1_sum_stats.append(metric1[0])
                g1_product_stats.append(metric1[1])
                g2_sum_stats.append(metric2[0])
                g2_product_stats.append(metric2[1])
            elif self.metric_type == "sum":
                g1_sum_stats.append(metric1)
                g2_sum_stats.append(metric2)
            else:
                g1_product_stats.append(metric1)
                g2_product_stats.append(metric2)
        if self.metric_type == "both":
            return (g1_sum_stats, g2_sum_stats), (g1_product_stats, g2_product_stats)
        elif self.metric_type == "sum":
            return g1_sum_stats, g2_sum_stats
        else:
            return g1_product_stats, g2_product_stats
            
    def differential_test(self, group1=None, group2=None, side="two-sided", seed=42):
        self._optional_print("Conducting differential analysis")

    
        if group1 is None or group2 is None:
            groups = self.data.iloc[:,0].unique()
            group_stats = []
            for group in groups:
                group_cells = self.data[self.data.iloc[:, 0] == group]
                num_cell = len(group_cells)
                if num_cell < self.MIN_GROUP_SIZE:
                    self._optional_print(f"Skipping group '{group}' (too few cells: {num_cell})")
                    continue
                    
                metric, num_cells = self._measure_one_group(group)
                if self.metric_type == "both":
                     stat_for_sort = metric[0]  
                else:
                     stat_for_sort = metric
                group_stats.append((group, stat_for_sort))
            group_stats = sorted(group_stats, key=lambda x: x[1], reverse=False)  
            group_pairs = list(zip([x[0] for x in group_stats[:-1]], [x[0] for x in group_stats[1:]]))
        else:
            if isinstance(group1, str): group1 = [group1]
            if isinstance(group2, str): group2 = [group2]
            if len(group1) != len(group2):
                raise ValueError("group1 and group2 should have the same sizes")
            group_pairs = list(zip(group1, group2))
    
        self._plot_groups = group_pairs

        out_df = pd.DataFrame(columns=[
        "groupA", "number of cells in groupA", "sum statistics for groupA", "product statistics for groupA",
        "groupB", "number of cells in groupB", "sum statistics for groupB","product statistics for groupB",
        "null hypothesis","test statistics sum", "test statistics product", "p value sum","p value product"
    ])
    
    
        for g1, g2 in group_pairs:
            p_sum, p_prod, null_hypothesis_str = self._pairwise_diff_test(g1, g2, seed=seed, side=side)

            metric_group1, group1_cells = self._measure_one_group(g1)
            metric_group2, group2_cells = self._measure_one_group(g2)

            if self.metric_type == "both":
                metric_sum_g1, metric_product_g1 = metric_group1
                metric_sum_g2, metric_product_g2 = metric_group2

                stat_sum = metric_sum_g1 - metric_sum_g2
                stat_product = metric_product_g1 - metric_product_g2

                new_row = {
                    "groupA": g1,
                    "number of cells in groupA": group1_cells,
                    "sum statistics for groupA": round(metric_sum_g1, DECIMAL),
                    "product statistics for groupA": round(metric_product_g1, DECIMAL),
                    "groupB": g2,
                    "number of cells in groupB": group2_cells,
                    "sum statistics for groupB": round(metric_sum_g2, DECIMAL),
                    "product statistics for groupB": round(metric_product_g2, DECIMAL),
                    "null hypothesis": null_hypothesis_str,
                    "test statistics sum": round(stat_sum, DECIMAL),
                    "test statistics product": round(stat_product, DECIMAL),
                    "p value sum": round(p_sum, DECIMAL),
                    "p value product": round(p_prod, DECIMAL)}
            else:
                stat = metric_group1 - metric_group2

                new_row = {
                    "groupA": g1,
                    "number of cells in groupA": group1_cells,
                    "sum statistics for groupA": round(metric_group1, DECIMAL) if self.metric_type=="sum" else None,
                    "product statistics for groupA": round(metric_group1, DECIMAL) if self.metric_type=="product" else None,
                    "groupB": g2,
                    "number of cells in groupB": group2_cells,
                    "sum statistics for groupB": round(metric_group2, DECIMAL) if self.metric_type=="sum" else None,
                    "product statistics for groupB": round(metric_group2, DECIMAL) if self.metric_type=="product" else None,
                    "null hypothesis": null_hypothesis_str,
                    "test statistics sum": round(stat, DECIMAL) if self.metric_type=="sum" else None,
                    "test statistics product": round(stat, DECIMAL) if self.metric_type=="product" else None,
                    "p value sum": round(p_sum, DECIMAL) if p_sum is not None else None,
                    "p value product": round(p_prod, DECIMAL) if p_prod is not None else None
                    
            }

            new_row_clean = {k:v for k,v in new_row.items() if pd.notna(v)}
            out_df = pd.concat([out_df, pd.DataFrame([new_row_clean])], ignore_index=True)
            out_df = out_df.dropna(axis=1, how='any')
    
        self._optional_print("Saving results")
        out_df.to_csv(os.path.join(self.output_dir, "differential_variability.csv"), index=False)
    
        return out_df
    
    def _pairwise_diff_test(self, group1, group2, a=DEFAULT_SIGNIFICANCE, seed=RANDOM_SEED, side="two-sided"):
        self._optional_print(f"Comparing {group1} and {group2}")
        random.seed(seed)
        self._dim_reduction()
        
        res = self._run_bootstrap(group1, group2, NUMBER_OF_BOOTSTRAP_SAMPLES, seed)
    
        if self.metric_type == "both":
            (g1_sum_stats, g2_sum_stats), (g1_product_stats, g2_product_stats) = res[0],res[1]

            p_val_sum = (np.sum(np.array(g1_sum_stats)[:, np.newaxis] >np.array(g2_sum_stats)) +1) / (len(np.array(g2_sum_stats)) * len(np.array(g1_sum_stats)) +1)
            p_val_product = (np.sum(np.array(g1_product_stats)[:, np.newaxis] >np.array(g2_product_stats)) +1) / (len(np.array(g2_product_stats)) * len(np.array(g1_product_stats)) +1)

            if side == "two-sided":
                p_val_sum = 2 * min(p_val_sum, 1-p_val_sum)
                p_val_product = 2 * min(p_val_product, 1-p_val_product)
                null_hypothesis_str = "A = B"
            elif side == "greater":
                p_val_sum = 1-p_val_sum
                p_val_product = 1-p_val_product
                null_hypothesis_str = "A < B"
            else:
                p_val_sum = p_val_sum
                p_val_product = p_val_product
                null_hypothesis_str = "A > B"
            return p_val_sum,p_val_product,null_hypothesis_str
                
        else:
            g1_stats, g2_stats = res
            diff = (np.array(g1_stats)[:,np.newaxis] - np.array(g2_stats)).flatten()
            p = (np.sum(np.array(g1_stats)[:,np.newaxis] >np.array(g2_stats)) +1) / (len(np.array(g1_stats)) * len(np.array(g2_stats)) +1)
            
            if side == "two-sided":
                p_val = 2 * min(p, 1-p)
                confidence_interval=(np.percentile(diff,100*0.5*0.05),np.percentile(diff,100-100*0.5*0.05))
                null_hypothesis_str = "A = B"
            elif side == "greater":
                p_val = 1-p
                confidence_interval = (np.percentile(diff,100*0.05), np.max(diff))
                null_hypothesis_str = "A < B"
            else:
                p_val = p
                confidence_interval = (np.min(diff), np.percentile(diff, 100*(1-0.05)))
                null_hypothesis_str = "A > B"
            return (p_val, None,null_hypothesis_str) if self.metric_type == "sum" else (None, p_val,null_hypothesis_str)


    def contain0(self,confidence_interval):
        min_value, max_value = confidence_interval
        return min_value < 0 < max_value
    
    
    def visualize_pca_results(self, detail=False, plot_pc_num=3,groups_to_plot=None):
        """visualization in pair-wiae PC space""" 
        if not detail:
            print("detail=False，skipping detailed PCA visualization.")
            return
        
        dimReduceData = self._dim_reduction(method="pca")

        print(self._plot_groups)
        if groups_to_plot is None:
            if hasattr(self, "_plot_groups"):
                group_pairs = self._plot_groups
            else:
                raise ValueError("No group pairs found. Please run differential_test first or pass group_pairs manually.")
        else:
            group_pairs=groups_to_plot
        
        selected_groups = set(itertools.chain(*group_pairs))  
        selected_groups = [g for g in selected_groups if g in dimReduceData.iloc[:, 0].unique()]
        plot_data = dimReduceData[dimReduceData.iloc[:, 0].isin(selected_groups)]
                
        summary_specific = []
        for group in selected_groups:
            group_data = plot_data[plot_data.iloc[:, 0] == group]
            group_matrix = self._change_of_axes(group_data)
            product_specific = self._summarize_cov_matrix(group_matrix.iloc[:, 1:])
            stds_specific = np.sqrt(np.diag(np.cov(group_matrix.iloc[:, 1:].T)))

            row_specific = [group]

            if self.metric_type in ["sum", "product"]:
                row_specific += stds_specific.tolist()[:self.num_pc]
                row_specific.append(product_specific)
            elif self.metric_type == "both":
                row_specific += stds_specific.tolist()[:self.num_pc]
                row_specific.append(product_specific)

            summary_specific.append(row_specific)



        col_names = ["Groups"]
        if self.metric_type == "sum":
            col_names += [f"PC{i+1}" for i in range(self.num_pc)]
            col_names += ["sum"]
        elif self.metric_type == "product":
            col_names += [f"PC{i+1}" for i in range(self.num_pc)]
            col_names += ["product"]
        elif self.metric_type == "both":
            col_names += [f"PC{i+1}" for i in range(self.num_pc)]
            col_names += ["sum","product"]

        summary_specific_df = pd.DataFrame(summary_specific, columns=col_names)

        if self.output_dir is not None:
            os.makedirs(self.output_dir, exist_ok=True)
            summary_specific_df.to_csv(os.path.join(self.output_dir, "summary_specific.csv"), index=False)

            
            
            
        pc_indices = list(range(1, plot_pc_num + 1))
        pc_pairs = list(itertools.combinations(pc_indices, 2))
        
        n_rows = len(group_pairs)
        n_cols = len(pc_pairs)
        fig_scat, axes_scat = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
        axes_scat = np.atleast_2d(axes_scat)
        fig_kde, axes_kde = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
        axes_kde = np.atleast_2d(axes_kde)
        
            
        for row_idx, (g1, g2) in enumerate(group_pairs):
            raw_data = self.data[self.data.iloc[:, 0].isin([g1, g2])]
            features = raw_data.iloc[:, 1:]
            n_components = min(self.num_pc, self.data.shape[1] - 1)
            pca = PCA(n_components)
            pcs = pca.fit_transform(features)
    
            plot_data = pd.concat([raw_data.iloc[:, 0].reset_index(drop=True),pd.DataFrame(pcs[:, :plot_pc_num], 
                            columns=[f"PC{i}" for i in range(1, plot_pc_num+1)])], axis=1)
    

            row_data = plot_data[[f"PC{i}" for i in range(1, plot_pc_num+1)]]
            min_val = row_data.min().min()
            max_val = row_data.max().max()
            margin = (max_val - min_val) * 0.1
            min_val -= margin
            max_val += margin

            groups = [g1, g2]
            colors = ["blue", "orange"]

            for col_idx, (pcx, pcy) in enumerate(pc_pairs):
                ax = axes_scat[row_idx, col_idx] if n_rows > 1 else axes_scat[col_idx]
                for g, c in zip(groups, colors):
                    subset = plot_data[plot_data.iloc[:, 0] == g]
                    ax.scatter(subset[f"PC{pcx}"], subset[f"PC{pcy}"], color=c, s=1, alpha=0.6, label=g)
                ax.set_xlabel(f"PC{pcx}")
                ax.set_ylabel(f"PC{pcy}")
                ax.set_xlim(min_val, max_val)
                ax.set_ylim(min_val, max_val)
                ax.set_aspect('equal', adjustable='box')
                if col_idx == 0:
                    ax.set_title(f"{g1} vs {g2}", fontsize=10)
                else:
                    ax.set_title("")
                ax.legend(groups, fontsize=6, loc='upper right', frameon=False)
            
            for col_idx, (pcx, pcy) in enumerate(pc_pairs):
                ax = axes_kde[row_idx, col_idx] if n_rows > 1 else axes_kde[col_idx]
                for g, c in zip(groups, colors):
                    subset = plot_data[plot_data.iloc[:, 0] == g]
                    if len(subset) > 1:
                        sns.kdeplot(
                    data=subset,
                    x=f"PC{pcx}", y=f"PC{pcy}",
                    fill=False, levels=8, linewidths=1, ax=ax, color=c
                    )
                ax.set_xlabel(f"PC{pcx}")
                ax.set_ylabel(f"PC{pcy}")
                ax.set_xlim(min_val, max_val)
                ax.set_ylim(min_val, max_val)
                ax.set_aspect('equal', adjustable='box')
                if col_idx == 0:
                    ax.set_title(f"{g1} vs {g2}", fontsize=10)
                else:
                    ax.set_title("")
                ax.legend(groups, fontsize=6, loc='upper right', frameon=False)

               

    
        if self.output_dir is not None:
            scatter_path = os.path.join(self.output_dir, "all_pairwise_scatter.png")
            kde_path = os.path.join(self.output_dir, "all_pairwise_kde.png")
            
            fig_scat.tight_layout(h_pad=2.0, w_pad=0.5)
            fig_scat.savefig(scatter_path, dpi=300)
            
            fig_kde.tight_layout(h_pad=2.0, w_pad=0.5)
            fig_kde.savefig(kde_path, dpi=300)

        plt.close(fig_scat)
        plt.close(fig_kde)
           
        print(f"All pairwise plots saved to {self.output_dir}/pairwise_kde_all.png" if self.output_dir else "Visualization complete.")
        
        
      
        