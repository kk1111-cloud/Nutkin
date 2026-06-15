import os
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import itertools
from .config import *
from sklearn.decomposition import PCA
import seaborn as sns
from scipy.stats import gaussian_kde

class Nutkin:
    def __init__(self, adata, output_path, num_pc = DEFAULT_NUM_PC, group_col="group", 
                 metric_type="sum", variability_method="mad", verbal=True, detail=False):
        self.verbal = verbal
        self.detail = detail
        self.group_col = group_col
        self.output_path = output_path
        self.num_pc = num_pc
        self.metric_type = metric_type.lower()
        self.variability_method = variability_method.lower()
        self.MIN_GROUP_SIZE = 3
                     
        if self.metric_type not in ["sum", "product", "both"]:
            raise ValueError("metric_type must be 'sum', 'product', or 'both'")
        
        if self.variability_method not in ["sd", "mad", "both"]:
            raise ValueError("variability_method must be 'sd' (standard deviation), 'mad' (median absolute deviation), or 'both'")
            
        os.makedirs(self.output_path, exist_ok=True)
        if group_col not in adata.obs.columns:
            raise ValueError(f"'{group_col}' not found in adata.obs.")
        self.adata = adata
        self.data = self.adata_to_df_with_obs(adata, group_col)
        self.dimReduceData = pd.DataFrame()
        
    def adata_to_df_with_obs(self, adata, group_col):
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
    
    def _optional_print(self, msg, *msgs):
        if self.verbal:
            additional_msg = ""
            for m in msgs:
                additional_msg += str(m) 
            combined_msg = msg + additional_msg
            print(combined_msg + '\n')
        return

    def _dim_reduction(self, method="pca"):
        if not self.dimReduceData.empty:
            return self.dimReduceData
        
        self._optional_print("Running dimensional reduction (PCA)...")
        n_components = min(self.num_pc, self.data.shape[1] - 1,self.data.shape[0] - 1)
        self.transformer = PCA(n_components)
        
        dimReduceData = pd.DataFrame(self.transformer.fit_transform(self.data.iloc[:, 1:]))
        self.dimReduceData = pd.concat([self.data.iloc[:,0], dimReduceData], axis=1)
        return self.dimReduceData

   
    def _change_of_axes(self, group_data):
        """represent group of cell into group specific features using pca"""
        n_components = min(self.num_pc, group_data.shape[1] - 1,group_data.shape[0] - 1)
        transformed_group_data = pd.DataFrame(PCA(n_components).fit_transform(group_data.iloc[:, 1:]))
        return pd.concat([group_data.iloc[:,0].reset_index(drop=True), transformed_group_data.reset_index(drop=True)], axis=1).reset_index(drop=True)
    
    
    def _summarize_cov_matrix(self, group_matrix):
        """Calculate variability metrics based on metric_type and variability_method"""
        eps = 1e-50
        
        # Calculate both methods if variability_method is "both"
        if self.variability_method in ["sd", "both"]:
            # Standard deviation method
            var_matrix = np.cov(group_matrix.T)
            
            if isinstance(var_matrix, float):
                diag_sd = np.array([var_matrix])
            else:
                diag_sd = np.diagonal(var_matrix).copy()
            
            diag_sd[diag_sd < eps] = eps
            
            metric_sum_sd = np.log(np.sum(np.sqrt(diag_sd)))
            metric_product_sd = np.sum(np.log(np.sqrt(diag_sd)))
        
        if self.variability_method in ["mad", "both"]:
            # Median absolute deviation method
            median_vec = np.median(group_matrix, axis=0)
            abs_dev = np.abs(group_matrix - median_vec)
            diag_mad = np.median(abs_dev, axis=0)
            
            diag_mad[diag_mad < eps] = eps
            
            metric_sum_mad = np.log(np.sum(diag_mad))
            metric_product_mad = np.sum(np.log(diag_mad))
        
        # Return based on variability_method and metric_type
        if self.variability_method == "sd":
            if self.metric_type == "sum":
                return metric_sum_sd
            elif self.metric_type == "product":
                return metric_product_sd
            else:  # metric_type == "both"
                return metric_sum_sd, metric_product_sd
                
        elif self.variability_method == "mad":
            if self.metric_type == "sum":
                return metric_sum_mad
            elif self.metric_type == "product":
                return metric_product_mad
            else:  # metric_type == "both"
                return metric_sum_mad, metric_product_mad
                
        else:  # variability_method == "both"
            if self.metric_type == "sum":
                return metric_sum_sd, metric_sum_mad
            elif self.metric_type == "product":
                return metric_product_sd, metric_product_mad
            else:  # metric_type == "both"
                return (metric_sum_sd, metric_product_sd), (metric_sum_mad, metric_product_mad)
    
    
    def measure(self):
        self._optional_print(f"Quantifying variability using {self.variability_method.upper()} method")
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
            
            # Build DataFrame based on variability_method and metric_type combinations
            if self.variability_method == "both":
                if self.metric_type == "both":
                    # 4 columns: sum_sd, product_sd, sum_mad, product_mad
                    (metric_sum_sd, metric_product_sd), (metric_sum_mad, metric_product_mad) = metric
                    sort_col = "sum statistics (sd)"
                    metric_df = pd.DataFrame([{
                        "group": group, 
                        "number of cells": num_cell, 
                        "sum statistics (sd)": round(metric_sum_sd, DECIMAL),
                        "product statistics (sd)": round(metric_product_sd, DECIMAL),
                        "sum statistics (mad)": round(metric_sum_mad, DECIMAL),
                        "product statistics (mad)": round(metric_product_mad, DECIMAL)
                    }])
                elif self.metric_type == "sum":
                    # 2 columns: sum_sd, sum_mad
                    metric_sum_sd, metric_sum_mad = metric
                    sort_col = "sum statistics (sd)"
                    metric_df = pd.DataFrame([{
                        "group": group,
                        "number of cells": num_cell,
                        "sum statistics (sd)": round(metric_sum_sd, DECIMAL),
                        "sum statistics (mad)": round(metric_sum_mad, DECIMAL)
                    }])
                else:  # product
                    # 2 columns: product_sd, product_mad
                    metric_product_sd, metric_product_mad = metric
                    sort_col = "product statistics (sd)"
                    metric_df = pd.DataFrame([{
                        "group": group,
                        "number of cells": num_cell,
                        "product statistics (sd)": round(metric_product_sd, DECIMAL),
                        "product statistics (mad)": round(metric_product_mad, DECIMAL)
                    }])
            elif self.metric_type == "both":
                # Single variability method, both metric types
                metric_sum, metric_product = metric
                sort_col = "sum statistics"
                metric_df = pd.DataFrame([{
                    "group": group, 
                    "number of cells": num_cell, 
                    "sum statistics": round(metric_sum, DECIMAL),
                    "product statistics": round(metric_product, DECIMAL)
                }])
            else:
                # Single variability method, single metric type
                sort_col = f"{self.metric_type} statistics"
                metric_df = pd.DataFrame([{
                    "group": group,
                    "number of cells": num_cell,
                    f"{self.metric_type} statistics": round(metric, DECIMAL)
                }])
                
            output_df = pd.concat([output_df, metric_df], ignore_index=True)
            
        output_df = output_df.sort_values(by=sort_col, ascending=True).reset_index(drop=True)
 
        self._optional_print("Saving results")
        output_df.to_csv(os.path.join(self.output_path, "var_measurement.csv"), index=False)
        return output_df
    
    def _measure_one_group(self, group):
        if self.dimReduceData.empty == True:
            self._dim_reduction()
        group_data = self.dimReduceData.loc[self.dimReduceData[self.dimReduceData.columns[0]] == group]
        group_matrix = self._change_of_axes(group_data)
        group_size = group_data.shape[0]
        metric = self._summarize_cov_matrix(group_matrix.iloc[:, 1:])
        return metric, group_size

    def _run_bootstrap(self, group1, group2, num_resample, seed=42):
        self._optional_print("Running Bootstrapping")
        random.seed(seed)
        self._dim_reduction()
        label_col_name = self.dimReduceData.columns.tolist()[0]
        group1_data = self.dimReduceData.loc[self.dimReduceData[label_col_name]==group1]
        group2_data = self.dimReduceData.loc[self.dimReduceData[label_col_name]==group2]
        
        num_resample_group2 = num_resample
        num_resample_group1 = num_resample
        
        # Initialize lists based on combinations
        if self.variability_method == "both":
            if self.metric_type == "both":
                g1_sum_sd, g2_sum_sd = [], []
                g1_product_sd, g2_product_sd = [], []
                g1_sum_mad, g2_sum_mad = [], []
                g1_product_mad, g2_product_mad = [], []
            elif self.metric_type == "sum":
                g1_sum_sd, g2_sum_sd = [], []
                g1_sum_mad, g2_sum_mad = [], []
            else:  # product
                g1_product_sd, g2_product_sd = [], []
                g1_product_mad, g2_product_mad = [], []
        elif self.metric_type == "both":
            g1_sum_stats, g2_sum_stats = [], []
            g1_product_stats, g2_product_stats = [], []
        else:
            g1_stats, g2_stats = [], []
        
        optimal_bootstrap_sample_size = min(BOOTSTRAP_SAMPLE_SIZE_UPPER_BOUND, int(MAXIMUM_BOOTSTRAP_PROPORTION*min(len(group1_data), len(group2_data))))
        
        for rep in range(num_resample):
            resample1 = group1_data.sample(optimal_bootstrap_sample_size, replace=True, random_state=seed + rep)
            resample2 = group2_data.sample(optimal_bootstrap_sample_size, replace=True, random_state=seed + rep)

            changed1 = self._change_of_axes(resample1)
            changed2 = self._change_of_axes(resample2)

            metric1 = self._summarize_cov_matrix(changed1.iloc[:, 1:])
            metric2 = self._summarize_cov_matrix(changed2.iloc[:, 1:])

            # Append based on combinations
            if self.variability_method == "both":
                if self.metric_type == "both":
                    (m1_sum_sd, m1_prod_sd), (m1_sum_mad, m1_prod_mad) = metric1
                    (m2_sum_sd, m2_prod_sd), (m2_sum_mad, m2_prod_mad) = metric2
                    g1_sum_sd.append(m1_sum_sd)
                    g1_product_sd.append(m1_prod_sd)
                    g1_sum_mad.append(m1_sum_mad)
                    g1_product_mad.append(m1_prod_mad)
                    g2_sum_sd.append(m2_sum_sd)
                    g2_product_sd.append(m2_prod_sd)
                    g2_sum_mad.append(m2_sum_mad)
                    g2_product_mad.append(m2_prod_mad)
                elif self.metric_type == "sum":
                    m1_sum_sd, m1_sum_mad = metric1
                    m2_sum_sd, m2_sum_mad = metric2
                    g1_sum_sd.append(m1_sum_sd)
                    g1_sum_mad.append(m1_sum_mad)
                    g2_sum_sd.append(m2_sum_sd)
                    g2_sum_mad.append(m2_sum_mad)
                else:  # product
                    m1_prod_sd, m1_prod_mad = metric1
                    m2_prod_sd, m2_prod_mad = metric2
                    g1_product_sd.append(m1_prod_sd)
                    g1_product_mad.append(m1_prod_mad)
                    g2_product_sd.append(m2_prod_sd)
                    g2_product_mad.append(m2_prod_mad)
            elif self.metric_type == "both":
                g1_sum_stats.append(metric1[0])
                g1_product_stats.append(metric1[1])
                g2_sum_stats.append(metric2[0])
                g2_product_stats.append(metric2[1])
            else:
                g1_stats.append(metric1)
                g2_stats.append(metric2)
        
        # Return based on combinations
        if self.variability_method == "both":
            if self.metric_type == "both":
                return ((g1_sum_sd, g2_sum_sd), (g1_product_sd, g2_product_sd)), ((g1_sum_mad, g2_sum_mad), (g1_product_mad, g2_product_mad))
            elif self.metric_type == "sum":
                return (g1_sum_sd, g2_sum_sd), (g1_sum_mad, g2_sum_mad)
            else:  # product
                return (g1_product_sd, g2_product_sd), (g1_product_mad, g2_product_mad)
        elif self.metric_type == "both":
            return (g1_sum_stats, g2_sum_stats), (g1_product_stats, g2_product_stats)
        else:
            return g1_stats, g2_stats
            
    def differential_test(self, group1=None, group2=None, side="two-sided", seed=42):
        self._optional_print(f"Conducting differential analysis using {self.variability_method.upper()} method")

    
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
                # Extract first metric for sorting
                if self.variability_method == "both":
                    if self.metric_type == "both":
                        stat_for_sort = metric[0][0]  # sum_sd
                    else:
                        stat_for_sort = metric[0]  # sd value
                elif self.metric_type == "both":
                    stat_for_sort = metric[0]  # sum
                else:
                    stat_for_sort = metric
                group_stats.append((group, stat_for_sort))
            group_stats = sorted(group_stats, key=lambda x: x[1], reverse=False)  
            group_pairs = list(zip([x[0] for x in group_stats[:-1]], [x[0] for x in group_stats[1:]]))
        else:
            if isinstance(group1, str): group1 = [group1]
            if isinstance(group2, str): group2 = [group2]
            print(group1)
            print(group2)
            if len(group1) != len(group2):
                raise ValueError("group1 and group2 should have the same sizes")
            group_pairs = list(zip(group1, group2))
    
        self._plot_groups = group_pairs

        # Build column names based on variability_method
        if self.variability_method == "both":
            out_df = pd.DataFrame(columns=[
                "groupA", "number of cells in groupA", 
                "sum statistics for groupA (sd)", "product statistics for groupA (sd)",
                "sum statistics for groupA (mad)", "product statistics for groupA (mad)",
                "groupB", "number of cells in groupB", 
                "sum statistics for groupB (sd)", "product statistics for groupB (sd)",
                "sum statistics for groupB (mad)", "product statistics for groupB (mad)",
                "null hypothesis",
                "test statistics sum (sd)", "test statistics product (sd)",
                "test statistics sum (mad)", "test statistics product (mad)",
                "p value sum (sd)", "p value product (sd)",
                "p value sum (mad)", "p value product (mad)"
            ])
        else:
            out_df = pd.DataFrame(columns=[
                "groupA", "number of cells in groupA", "sum statistics for groupA", "product statistics for groupA",
                "groupB", "number of cells in groupB", "sum statistics for groupB","product statistics for groupB",
                "null hypothesis","test statistics sum", "test statistics product", "p value sum","p value product"
            ])
    
    
        for g1, g2 in group_pairs:
            p_results, null_hypothesis_str = self._pairwise_diff_test(g1, g2, seed=seed, side=side)

            metric_group1, group1_cells = self._measure_one_group(g1)
            metric_group2, group2_cells = self._measure_one_group(g2)

            if self.variability_method == "both":
                if self.metric_type == "both":
                    (m1_sum_sd, m1_prod_sd), (m1_sum_mad, m1_prod_mad) = metric_group1
                    (m2_sum_sd, m2_prod_sd), (m2_sum_mad, m2_prod_mad) = metric_group2
                    p_sum_sd, p_prod_sd, p_sum_mad, p_prod_mad = p_results
                    
                    new_row = {
                        "groupA": g1,
                        "number of cells in groupA": group1_cells,
                        "sum statistics for groupA (sd)": round(m1_sum_sd, DECIMAL),
                        "product statistics for groupA (sd)": round(m1_prod_sd, DECIMAL),
                        "sum statistics for groupA (mad)": round(m1_sum_mad, DECIMAL),
                        "product statistics for groupA (mad)": round(m1_prod_mad, DECIMAL),
                        "groupB": g2,
                        "number of cells in groupB": group2_cells,
                        "sum statistics for groupB (sd)": round(m2_sum_sd, DECIMAL),
                        "product statistics for groupB (sd)": round(m2_prod_sd, DECIMAL),
                        "sum statistics for groupB (mad)": round(m2_sum_mad, DECIMAL),
                        "product statistics for groupB (mad)": round(m2_prod_mad, DECIMAL),
                        "null hypothesis": null_hypothesis_str,
                        "test statistics sum (sd)": round(m1_sum_sd - m2_sum_sd, DECIMAL),
                        "test statistics product (sd)": round(m1_prod_sd - m2_prod_sd, DECIMAL),
                        "test statistics sum (mad)": round(m1_sum_mad - m2_sum_mad, DECIMAL),
                        "test statistics product (mad)": round(m1_prod_mad - m2_prod_mad, DECIMAL),
                        "p value sum (sd)": round(p_sum_sd, DECIMAL),
                        "p value product (sd)": round(p_prod_sd, DECIMAL),
                        "p value sum (mad)": round(p_sum_mad, DECIMAL),
                        "p value product (mad)": round(p_prod_mad, DECIMAL)
                    }
                elif self.metric_type == "sum":
                    m1_sum_sd, m1_sum_mad = metric_group1
                    m2_sum_sd, m2_sum_mad = metric_group2
                    p_sum_sd, p_sum_mad = p_results
                    
                    new_row = {
                        "groupA": g1,
                        "number of cells in groupA": group1_cells,
                        "sum statistics for groupA (sd)": round(m1_sum_sd, DECIMAL),
                        "sum statistics for groupA (mad)": round(m1_sum_mad, DECIMAL),
                        "groupB": g2,
                        "number of cells in groupB": group2_cells,
                        "sum statistics for groupB (sd)": round(m2_sum_sd, DECIMAL),
                        "sum statistics for groupB (mad)": round(m2_sum_mad, DECIMAL),
                        "null hypothesis": null_hypothesis_str,
                        "test statistics sum (sd)": round(m1_sum_sd - m2_sum_sd, DECIMAL),
                        "test statistics sum (mad)": round(m1_sum_mad - m2_sum_mad, DECIMAL),
                        "p value sum (sd)": round(p_sum_sd, DECIMAL),
                        "p value sum (mad)": round(p_sum_mad, DECIMAL)
                    }
                else:  # product
                    m1_prod_sd, m1_prod_mad = metric_group1
                    m2_prod_sd, m2_prod_mad = metric_group2
                    p_prod_sd, p_prod_mad = p_results
                    
                    new_row = {
                        "groupA": g1,
                        "number of cells in groupA": group1_cells,
                        "product statistics for groupA (sd)": round(m1_prod_sd, DECIMAL),
                        "product statistics for groupA (mad)": round(m1_prod_mad, DECIMAL),
                        "groupB": g2,
                        "number of cells in groupB": group2_cells,
                        "product statistics for groupB (sd)": round(m2_prod_sd, DECIMAL),
                        "product statistics for groupB (mad)": round(m2_prod_mad, DECIMAL),
                        "null hypothesis": null_hypothesis_str,
                        "test statistics product (sd)": round(m1_prod_sd - m2_prod_sd, DECIMAL),
                        "test statistics product (mad)": round(m1_prod_mad - m2_prod_mad, DECIMAL),
                        "p value product (sd)": round(p_prod_sd, DECIMAL),
                        "p value product (mad)": round(p_prod_mad, DECIMAL)
                    }
            elif self.metric_type == "both":
                metric_sum_g1, metric_product_g1 = metric_group1
                metric_sum_g2, metric_product_g2 = metric_group2
                p_sum, p_prod = p_results

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
                p_val = p_results[0] if isinstance(p_results, tuple) else p_results

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
                    "p value sum": round(p_val, DECIMAL) if self.metric_type=="sum" else None,
                    "p value product": round(p_val, DECIMAL) if self.metric_type=="product" else None
                    
            }

            new_row_clean = {k:v for k,v in new_row.items() if pd.notna(v)}
            out_df = pd.concat([out_df, pd.DataFrame([new_row_clean])], ignore_index=True)
            
        out_df = out_df.dropna(axis=1, how='all')
    
        self._optional_print("Saving results")
        out_df.to_csv(os.path.join(self.output_path, "differential_variability.csv"), index=False)
    
        return out_df
    
    def _pairwise_diff_test(self, group1, group2, a=DEFAULT_SIGNIFICANCE, seed=RANDOM_SEED, side="two-sided"):
        self._optional_print(f"Comparing {group1} and {group2}")
        random.seed(seed)
        self._dim_reduction()
        
        res = self._run_bootstrap(group1, group2, NUMBER_OF_BOOTSTRAP_SAMPLES, seed)
    
        if self.variability_method == "both":
            if self.metric_type == "both":
                ((g1_sum_sd, g2_sum_sd), (g1_prod_sd, g2_prod_sd)), ((g1_sum_mad, g2_sum_mad), (g1_prod_mad, g2_prod_mad)) = res
                
                p_sum_sd = (np.sum(np.array(g1_sum_sd)[:, np.newaxis] > np.array(g2_sum_sd)) + 1) / (len(g2_sum_sd) * len(g1_sum_sd) + 1)
                p_prod_sd = (np.sum(np.array(g1_prod_sd)[:, np.newaxis] > np.array(g2_prod_sd)) + 1) / (len(g2_prod_sd) * len(g1_prod_sd) + 1)
                p_sum_mad = (np.sum(np.array(g1_sum_mad)[:, np.newaxis] > np.array(g2_sum_mad)) + 1) / (len(g2_sum_mad) * len(g1_sum_mad) + 1)
                p_prod_mad = (np.sum(np.array(g1_prod_mad)[:, np.newaxis] > np.array(g2_prod_mad)) + 1) / (len(g2_prod_mad) * len(g1_prod_mad) + 1)
                
                if side == "two-sided":
                    p_sum_sd = 2 * min(p_sum_sd, 1 - p_sum_sd)
                    p_prod_sd = 2 * min(p_prod_sd, 1 - p_prod_sd)
                    p_sum_mad = 2 * min(p_sum_mad, 1 - p_sum_mad)
                    p_prod_mad = 2 * min(p_prod_mad, 1 - p_prod_mad)
                    null_hypothesis_str = "A = B"
                elif side == "greater":
                    p_sum_sd = 1 - p_sum_sd
                    p_prod_sd = 1 - p_prod_sd
                    p_sum_mad = 1 - p_sum_mad
                    p_prod_mad = 1 - p_prod_mad
                    null_hypothesis_str = "A < B"
                else:
                    null_hypothesis_str = "A > B"
                
                return (p_sum_sd, p_prod_sd, p_sum_mad, p_prod_mad), null_hypothesis_str
                
            elif self.metric_type == "sum":
                (g1_sum_sd, g2_sum_sd), (g1_sum_mad, g2_sum_mad) = res
                
                p_sum_sd = (np.sum(np.array(g1_sum_sd)[:, np.newaxis] > np.array(g2_sum_sd)) + 1) / (len(g2_sum_sd) * len(g1_sum_sd) + 1)
                p_sum_mad = (np.sum(np.array(g1_sum_mad)[:, np.newaxis] > np.array(g2_sum_mad)) + 1) / (len(g2_sum_mad) * len(g1_sum_mad) + 1)
                
                if side == "two-sided":
                    p_sum_sd = 2 * min(p_sum_sd, 1 - p_sum_sd)
                    p_sum_mad = 2 * min(p_sum_mad, 1 - p_sum_mad)
                    null_hypothesis_str = "A = B"
                elif side == "greater":
                    p_sum_sd = 1 - p_sum_sd
                    p_sum_mad = 1 - p_sum_mad
                    null_hypothesis_str = "A < B"
                else:
                    null_hypothesis_str = "A > B"
                
                return (p_sum_sd, p_sum_mad), null_hypothesis_str
                
            else:  # product
                (g1_prod_sd, g2_prod_sd), (g1_prod_mad, g2_prod_mad) = res
                
                p_prod_sd = (np.sum(np.array(g1_prod_sd)[:, np.newaxis] > np.array(g2_prod_sd)) + 1) / (len(g2_prod_sd) * len(g1_prod_sd) + 1)
                p_prod_mad = (np.sum(np.array(g1_prod_mad)[:, np.newaxis] > np.array(g2_prod_mad)) + 1) / (len(g2_prod_mad) * len(g1_prod_mad) + 1)
                
                if side == "two-sided":
                    p_prod_sd = 2 * min(p_prod_sd, 1 - p_prod_sd)
                    p_prod_mad = 2 * min(p_prod_mad, 1 - p_prod_mad)
                    null_hypothesis_str = "A = B"
                elif side == "greater":
                    p_prod_sd = 1 - p_prod_sd
                    p_prod_mad = 1 - p_prod_mad
                    null_hypothesis_str = "A < B"
                else:
                    null_hypothesis_str = "A > B"
                
                return (p_prod_sd, p_prod_mad), null_hypothesis_str
                
        elif self.metric_type == "both":
            (g1_sum_stats, g2_sum_stats), (g1_product_stats, g2_product_stats) = res

            p_val_sum = (np.sum(np.array(g1_sum_stats)[:, np.newaxis] > np.array(g2_sum_stats)) + 1) / (len(np.array(g2_sum_stats)) * len(np.array(g1_sum_stats)) + 1)
            p_val_product = (np.sum(np.array(g1_product_stats)[:, np.newaxis] > np.array(g2_product_stats)) + 1) / (len(np.array(g2_product_stats)) * len(np.array(g1_product_stats)) + 1)

            if side == "two-sided":
                p_val_sum = 2 * min(p_val_sum, 1 - p_val_sum)
                p_val_product = 2 * min(p_val_product, 1 - p_val_product)
                null_hypothesis_str = "A = B"
            elif side == "greater":
                p_val_sum = 1 - p_val_sum
                p_val_product = 1 - p_val_product
                null_hypothesis_str = "A < B"
            else:
                null_hypothesis_str = "A > B"
            return (p_val_sum, p_val_product), null_hypothesis_str
                
        else:
            g1_stats, g2_stats = res
            diff = (np.array(g1_stats)[:, np.newaxis] - np.array(g2_stats)).flatten()
            p = (np.sum(np.array(g1_stats)[:, np.newaxis] > np.array(g2_stats)) + 1) / (len(np.array(g1_stats)) * len(np.array(g2_stats)) + 1)
            
            if side == "two-sided":
                p_val = 2 * min(p, 1 - p)
                confidence_interval = (np.percentile(diff, 100 * 0.5 * 0.05), np.percentile(diff, 100 - 100 * 0.5 * 0.05))
                null_hypothesis_str = "A = B"
            elif side == "greater":
                p_val = 1 - p
                confidence_interval = (np.percentile(diff, 100 * 0.05), np.max(diff))
                null_hypothesis_str = "A < B"
            else:
                p_val = p
                confidence_interval = (np.min(diff), np.percentile(diff, 100 * (1 - 0.05)))
                null_hypothesis_str = "A > B"
            return (p_val,), null_hypothesis_str


    def contain0(self, confidence_interval):
        min_value, max_value = confidence_interval
        return min_value < 0 < max_value
    
    
    def visualize_pca_results(self, detail=False, plot_pc_num=3, plot_type="scatter", groups_to_plot=None):
        """visualization in pair-wise PC space""" 
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
            metric_result = self._summarize_cov_matrix(group_matrix.iloc[:, 1:])
            stds_specific = np.sqrt(np.diag(np.cov(group_matrix.iloc[:, 1:].T)))

            row_specific = [group]
            row_specific += stds_specific.tolist()[:self.num_pc]
            
            # Handle different combinations of variability_method and metric_type
            if self.variability_method == "both":
                if self.metric_type == "both":
                    (sum_sd, prod_sd), (sum_mad, prod_mad) = metric_result
                    row_specific += [sum_sd, prod_sd, sum_mad, prod_mad]
                elif self.metric_type == "sum":
                    sum_sd, sum_mad = metric_result
                    row_specific += [sum_sd, sum_mad]
                else:  # product
                    prod_sd, prod_mad = metric_result
                    row_specific += [prod_sd, prod_mad]
            elif self.metric_type == "both":
                sum_stat, prod_stat = metric_result
                row_specific += [sum_stat, prod_stat]
            else:
                row_specific.append(metric_result)

            summary_specific.append(row_specific)

        # Build column names based on variability_method and metric_type
        col_names = ["Groups"] + [f"PC{i+1}" for i in range(self.num_pc)]
        
        if self.variability_method == "both":
            if self.metric_type == "both":
                col_names += ["sum (sd)", "product (sd)", "sum (mad)", "product (mad)"]
            elif self.metric_type == "sum":
                col_names += ["sum (sd)", "sum (mad)"]
            else:  # product
                col_names += ["product (sd)", "product (mad)"]
        elif self.metric_type == "both":
            col_names += ["sum", "product"]
        elif self.metric_type == "sum":
            col_names += ["sum"]
        else:  # product
            col_names += ["product"]

        summary_specific_df = pd.DataFrame(summary_specific, columns=col_names)

        if self.output_path is not None:
            os.makedirs(self.output_path, exist_ok=True)
            summary_specific_df.to_csv(os.path.join(self.output_path, "summary_specific.csv"), index=False)

        # Prepare data for plotting
        X = self.data.iloc[:, 1:].values
        groups_all = self.data.iloc[:, 0].values

        pca_global = PCA(n_components=plot_pc_num)
        pcs_global = pca_global.fit_transform(X)
        pca_df = pd.DataFrame(
            pcs_global,
            columns=[f"PC{i}" for i in range(1, plot_pc_num + 1)]
        )
        pca_df["group"] = groups_all
            
        plot_data = []
            
        for g in selected_groups:
            df_g = pca_df[pca_df["group"] == g].copy()
            if len(df_g) < 2:
                continue
            pca_g = PCA(n_components=plot_pc_num)
            pcs_g = pca_g.fit_transform(
                df_g[[f"PC{i}" for i in range(1, plot_pc_num + 1)]]
            )
            
            for i in range(plot_pc_num):
                df_g[f"subPC{i+1}"] = pcs_g[:, i]

            plot_data.append(df_g)

        plot_data = pd.concat(plot_data, ignore_index=True)

        # Generate all PC pairs
        pc_indices = list(range(1, plot_pc_num + 1))
        pc_pairs = list(itertools.combinations(pc_indices, 2))

        n_rows = len(group_pairs)
        n_cols = len(pc_pairs)

        # Normalize plot_type to list
        if isinstance(plot_type, str):
            plot_types = [plot_type]
        else:
            plot_types = plot_type

        # Validate plot types
        valid_types = {"scatter", "contour", "both"}
        for pt in plot_types:
            if pt not in valid_types:
                raise ValueError(f"plot_type must be 'scatter', 'contour', 'both', or a list of these. Got: {pt}")

        colors = ["blue", "orange"]

        # Create separate figures for each plot type
        for current_plot_type in plot_types:
            fig, axes = plt.subplots(
                n_rows, n_cols,
                figsize=(4 * n_cols, 3 * n_rows)
            )
            axes = np.atleast_2d(axes)

            for row_idx, (g1, g2) in enumerate(group_pairs):
                for col_idx, (pcx, pcy) in enumerate(pc_pairs):
                    ax = axes[row_idx, col_idx]

                    for g, c in zip([g1, g2], colors):
                        subset = plot_data[plot_data["group"] == g]

                        if current_plot_type == "scatter":
                            ax.scatter(
                                subset[f"subPC{pcx}"],
                                subset[f"subPC{pcy}"],
                                s=2,
                                alpha=0.6,
                                color=c,
                                label=g
                            )

                        elif current_plot_type == "contour":
                            self._draw_kde_contour(ax, subset, pcx, pcy, c)

                        elif current_plot_type == "both":
                            # First draw contour, then scatter on top
                            self._draw_kde_contour(ax, subset, pcx, pcy, c)
                            ax.scatter(
                                subset[f"subPC{pcx}"],
                                subset[f"subPC{pcy}"],
                                s=2,
                                alpha=0.4,
                                color=c,
                                label=g
                            )

                    ax.set_xlabel(f"subPC{pcx}")
                    ax.set_ylabel(f"subPC{pcy}")
                    ax.set_aspect("equal", adjustable="box")

                    if col_idx == 0:
                        ax.set_title(f"{g1} vs {g2}", fontsize=10)

                    ax.legend(fontsize=6, frameon=False)

            plt.tight_layout()

            if self.output_path is not None:
                fname = f"pairwise_pca_{current_plot_type}.png"
                fig.savefig(
                    os.path.join(self.output_path, fname),
                    dpi=300,
                    bbox_inches="tight"
                )

            plt.close(fig)

            print(
                f"PCA visualization ({current_plot_type}) saved to {self.output_path}"
                if self.output_path else
                f"Visualization ({current_plot_type}) complete."
            )
    
    def _draw_kde_contour(self, ax, subset, pcx, pcy, color):
        """Helper function to draw KDE contour plots"""
        x = subset[f"subPC{pcx}"].values
        y = subset[f"subPC{pcy}"].values

        if len(x) < 10:
            return

        kde = gaussian_kde(np.vstack([x, y]))

        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        pad_x = 0.1 * (xmax - xmin) if xmax != xmin else 0.1
        pad_y = 0.1 * (ymax - ymin) if ymax != ymin else 0.1

        xx, yy = np.mgrid[
            xmin - pad_x : xmax + pad_x : 100j,
            ymin - pad_y : ymax + pad_y : 100j
        ]

        zz = kde(
            np.vstack([xx.ravel(), yy.ravel()])
        ).reshape(xx.shape)

        # Low-density shape
        ax.contour(
            xx, yy, zz,
            levels=10,
            colors=color,
            linewidths=0.8,
            alpha=0.8
        )

        # Fixed density levels (if defined in config)
        try:
            from config import contour_levels_fixed
            cs = ax.contour(
                xx, yy, zz,
                levels=contour_levels_fixed,
                colors=color,
                linewidths=1.6
            )
            ax.clabel(cs, fmt="%.1e", fontsize=8)
        except (ImportError, AttributeError):
            # If contour_levels_fixed not defined, skip this part
            pass
