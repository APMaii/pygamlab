import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr, kendalltau, jarque_bera, normaltest
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.impute import KNNImputer
from datetime import datetime
import warnings
#warnings.filterwarnings('ignore')
#plt.style.use('default')
#sns.set_palette("viridis")
#sns.set_context("talk")



class Automated_Exploratory_Processor:
    """
    A comprehensive automated exploratory data analysis tool for any dataset
    Automatically detects numeric columns and performs statistical analysis
    """
    
    def __init__(self, data_source, numeric_columns=None):
        """
        Initialize the processor with data
        
        Parameters:
        -----------
        data_source : str or pd.DataFrame
            Path to Excel file (str) or pandas DataFrame
        numeric_columns : list, optional
            List of column names to analyze. If None, automatically detects numeric columns
        """
        # Load data from Excel file or use DataFrame directly
        if isinstance(data_source, str):
            if 'xlsx' in data_source or 'xls' in data_source:
                self.df = pd.read_excel(data_source)
            elif 'csv' in data_source:
                self.df = pd.read_csv(data_source)
            elif 'txt' in data_source:
                self.df = pd.read_csv(data_source, sep='\t')
            elif 'json' in data_source:
                self.df = pd.read_json(data_source)
            else:
                raise ValueError("Unsupported file format. Please provide a valid Excel file.")
        elif isinstance(data_source, pd.DataFrame):
            self.df = data_source.copy()
        else:
            raise ValueError("data_source must be a file path (str) or pandas DataFrame")
        
        # Auto-detect numeric columns if not specified
        if numeric_columns is None:
            self.numeric_columns = self.df.select_dtypes(include=[np.number]).columns.tolist()
        else:
            # Validate that specified columns exist and are numeric
            available_cols = [c for c in numeric_columns if c in self.df.columns]
            numeric_cols = [c for c in available_cols if pd.api.types.is_numeric_dtype(self.df[c])]
            self.numeric_columns = numeric_cols
        
        if len(self.numeric_columns) == 0:
            raise ValueError("No numeric columns found in the dataset")
        
        print(f"Dataset loaded: {self.df.shape[0]} rows, {self.df.shape[1]} columns")
        print(f"Numeric columns detected: {self.numeric_columns}")


    def preprocess_outliers(self, target_column=None, method='double_iqr'):
        """
        Remove outliers from a target column using IQR method
        
        Parameters:
        -----------
        target_column : str, optional
            Column name to clean. If None, uses first numeric column
        method : str
            'single_iqr' or 'double_iqr' (default: double IQR cleaning)
        
        Returns:
        --------
        pd.DataFrame, int : cleaned dataframe and final length
        """
        if target_column is None:
            target_column = self.numeric_columns[0]
        
        print("\n" + "="*80)
        print(f"PREPROCESSING: Outlier Removal for '{target_column}'")
        print("="*80)

        if target_column not in self.df.columns:
            print(f"Column '{target_column}' not found. Skipping preprocessing.")
            return self.df, len(self.df)
        
        if target_column not in self.numeric_columns:
            print(f"Column '{target_column}' is not numeric. Skipping preprocessing.")
            return self.df, len(self.df)
        
        initial_len = len(self.df)

        # First IQR cleaning
        Q1 = np.percentile(self.df[target_column].dropna(), 25)
        Q3 = np.percentile(self.df[target_column].dropna(), 75)
        IQR = Q3 - Q1
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR
        
        mask1 = (self.df[target_column] >= lower) & (self.df[target_column] <= upper)
        df_clean = self.df[mask1].copy()
        
        print(f"Step 1 removed: {initial_len - len(df_clean)} outliers")

        if method == 'double_iqr':
            # Second IQR cleaning on the cleaned set
            Q1_2 = np.percentile(df_clean[target_column].dropna(), 25)
            Q3_2 = np.percentile(df_clean[target_column].dropna(), 75)
            IQR_2 = Q3_2 - Q1_2
            lower_2 = Q1_2 - 1.5 * IQR_2
            upper_2 = Q3_2 + 1.5 * IQR_2
            
            mask2 = (df_clean[target_column] >= lower_2) & (df_clean[target_column] <= upper_2)
            df_clean2 = df_clean[mask2].copy()
            
            print(f"Step 2 removed: {len(df_clean) - len(df_clean2)} outliers")
            total_removed = initial_len - len(df_clean2)
            final_df = df_clean2
        else:
            total_removed = initial_len - len(df_clean)
            final_df = df_clean
        
        print("-"*80)
        print(f"Total removed outliers: {total_removed}")
        print(f"Final dataset length: {len(final_df)}")

        # Update original dataframe
        self.df = final_df.reset_index(drop=True)
        return self.df, len(self.df)


    def basic_statistics(self):
        print("\n" + "="*70)
        print("BASIC STATISTICS")
        print("="*70)
        print(f"Shape: {self.df.shape}")
        print("\nMissing values:")
        print(self.df[self.numeric_columns].isnull().sum())
        stats_df = self.df[self.numeric_columns].describe()
        print(stats_df)
        return stats_df
    

    def correlation_analysis(self):
        print("\n" + "="*70)
        print("CORRELATION MATRIX")
        print("="*70)
        pearson_corr = self.df[self.numeric_columns].corr()
        print(pearson_corr.round(3))

        strong_corrs = []
        for i in range(len(pearson_corr.columns)):
            for j in range(i+1, len(pearson_corr.columns)):
                if abs(pearson_corr.iloc[i,j]) > 0.7:
                    strong_corrs.append((pearson_corr.columns[i],
                                         pearson_corr.columns[j],
                                         pearson_corr.iloc[i,j]))

        print("\nStrong correlations (|r| > 0.7):")
        for a,b,r in strong_corrs:
            print(f"{a} vs {b}: r = {r:.3f}")

        return pearson_corr
    

    def plot_correlation_heatmap(self):
        corr = self.df[self.numeric_columns].corr()
        plt.figure(figsize=(10, 8))
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, annot=True, cmap="coolwarm", vmin=-1, vmax=1,
                    linewidths=1, square=True, mask=mask)
        plt.title("Correlation Heatmap", fontsize=18)
        plt.tight_layout()
        plt.savefig("correlation_heatmap.png", dpi=300)
        print("✓ Saved: correlation_heatmap.png")
        plt.show()
    

    def plot_distributions(self):
        n = len(self.numeric_columns)
        rows = (n + 2) // 3
        plt.figure(figsize=(18, 5*rows))

        for i, col in enumerate(self.numeric_columns, 1):
            plt.subplot(rows, 3, i)
            sns.histplot(self.df[col], kde=True, edgecolor="black")
            plt.title(f"{col} Distribution")
            plt.xlabel(col)
        plt.tight_layout()
        plt.savefig("distributions.png", dpi=300)
        print("✓ Saved: distributions.png")
        plt.show()


    def normality_tests(self):
        print("\n" + "="*80)
        print("NORMALITY TESTS (Shapiro-Wilk)")
        print("="*80)
        print("H0: Data is normally distributed (p > 0.05)")
        print("-"*80)
        
        for col in self.numeric_columns:
            stat, p_value = stats.shapiro(self.df[col].dropna())
            result = "Normal" if p_value > 0.05 else "Not Normal"
            print(f"{col:20s}: statistic = {stat:.4f}, p-value = {p_value:.4f} [{result}]")
    

    def outlier_detection(self):
        print("\n" + "="*80)
        print("OUTLIER DETECTION")
        print("="*80)
        
        outliers_summary = {}
        
        for col in self.numeric_columns:
            Q1 = self.df[col].quantile(0.25)
            Q3 = self.df[col].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            
            iqr_outliers = self.df[(self.df[col] < lower_bound) | 
                                   (self.df[col] > upper_bound)]
            
            z_scores = np.abs(stats.zscore(self.df[col].dropna()))
            z_outliers = self.df[col][z_scores > 3]
            
            outliers_summary[col] = {
                'IQR_outliers': len(iqr_outliers),
                'Z_score_outliers': len(z_outliers)
            }
            
            if len(iqr_outliers) > 0:
                print(f"\n{col}:")
                print(f"  IQR outliers: {len(iqr_outliers)} samples")
                print(f"  Z-score outliers: {len(z_outliers)} samples")
        
        return outliers_summary
    

    def boxplots(self):
        n = len(self.numeric_columns)
        rows = (n + 2) // 3
        plt.figure(figsize=(18, 5*rows))

        for i, col in enumerate(self.numeric_columns, 1):
            plt.subplot(rows, 3, i)
            sns.boxplot(y=self.df[col])
            plt.title(f"{col}")
        plt.tight_layout()
        plt.savefig("boxplots.png", dpi=300)
        print("✓ Saved: boxplots.png")
        plt.show()
    

    def plot_pairwise_relationships(self, target=None):
        """
        Plot pairwise relationships between target and other numeric variables
        
        Parameters:
        -----------
        target : str, optional
            Target column name. If None, uses first numeric column
        """
        if target is None:
            target = self.numeric_columns[0]
        
        if target not in self.numeric_columns:
            print(f"Target '{target}' not found in numeric columns.")
            return

        vars_other = [c for c in self.numeric_columns if c != target]
        if len(vars_other) == 0:
            print("No other numeric columns to compare with target.")
            return
        
        n = len(vars_other)
        rows = (n + 2) // 3
        plt.figure(figsize=(18, 5*rows))

        for i, var in enumerate(vars_other, 1):
            plt.subplot(rows, 3, i)
            sns.regplot(x=self.df[var], y=self.df[target],
                        scatter_kws={'s':50, 'edgecolor':'black'})
            corr, _ = pearsonr(self.df[var].dropna(), self.df[target].dropna())
            plt.title(f"{var} vs {target}\nr={corr:.2f}")
        
        plt.tight_layout()
        plt.savefig(f"{target}_relationships.png", dpi=300)
        print(f"✓ Saved: {target}_relationships.png")
        plt.show()
    

    def missing_data_analysis(self):
        """
        Comprehensive missing data analysis
        """
        print("\n" + "="*80)
        print("MISSING DATA ANALYSIS")
        print("="*80)
        
        missing_count = self.df[self.numeric_columns].isnull().sum()
        missing_pct = (missing_count / len(self.df)) * 100
        
        missing_df = pd.DataFrame({
            'Column': missing_count.index,
            'Missing_Count': missing_count.values,
            'Missing_Percentage': missing_pct.values
        })
        missing_df = missing_df[missing_df['Missing_Count'] > 0].sort_values('Missing_Percentage', ascending=False)
        
        if len(missing_df) > 0:
            print("\nMissing Data Summary:")
            print(missing_df.to_string(index=False))
            
            # Visualize missing data
            plt.figure(figsize=(12, 6))
            sns.barplot(data=missing_df, x='Column', y='Missing_Percentage', palette='viridis')
            plt.title('Missing Data Percentage by Column', fontsize=16)
            plt.xlabel('Column')
            plt.ylabel('Missing Percentage (%)')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig("missing_data_analysis.png", dpi=300)
            print("✓ Saved: missing_data_analysis.png")
            plt.show()
        else:
            print("✓ No missing values found in numeric columns")
        
        return missing_df
    
    
    def impute_missing_values(self, strategy='mean', columns=None):
        """
        Impute missing values using various strategies
        
        Parameters:
        -----------
        strategy : str
            'mean', 'median', 'mode', 'knn', 'forward_fill', 'backward_fill'
        columns : list, optional
            Columns to impute. If None, imputes all numeric columns
        """
        if columns is None:
            columns = self.numeric_columns
        
        print(f"\nImputing missing values using '{strategy}' strategy...")
        
        df_imputed = self.df.copy()
        
        for col in columns:
            if col not in self.df.columns:
                continue
            
            missing_count = self.df[col].isnull().sum()
            if missing_count == 0:
                continue
            
            if strategy == 'mean':
                df_imputed[col].fillna(self.df[col].mean(), inplace=True)
            elif strategy == 'median':
                df_imputed[col].fillna(self.df[col].median(), inplace=True)
            elif strategy == 'mode':
                df_imputed[col].fillna(self.df[col].mode()[0], inplace=True)
            elif strategy == 'forward_fill':
                df_imputed[col] = df_imputed[col].ffill()
            elif strategy == 'backward_fill':
                df_imputed[col] = df_imputed[col].bfill()
            elif strategy == 'knn':
                imputer = KNNImputer(n_neighbors=5)
                df_imputed[[col]] = imputer.fit_transform(self.df[[col]])
        
        self.df = df_imputed
        print(f"✓ Imputation complete. Remaining missing values: {self.df[columns].isnull().sum().sum()}")
        return self.df
    
    
    def skewness_kurtosis_analysis(self):
        """
        Analyze skewness and kurtosis for all numeric columns
        """
        print("\n" + "="*80)
        print("SKEWNESS & KURTOSIS ANALYSIS")
        print("="*80)
        print("Skewness: |skew| < 0.5 (approximately symmetric), |skew| > 1 (highly skewed)")
        print("Kurtosis: ~3 (normal), >3 (heavy tails), <3 (light tails)")
        print("-"*80)
        
        results = []
        for col in self.numeric_columns:
            skew = stats.skew(self.df[col].dropna())
            kurt = stats.kurtosis(self.df[col].dropna())
            
            skew_interpretation = "Symmetric" if abs(skew) < 0.5 else ("Moderately Skewed" if abs(skew) < 1 else "Highly Skewed")
            kurt_interpretation = "Normal" if abs(kurt - 3) < 0.5 else ("Heavy Tails" if kurt > 3 else "Light Tails")
            
            results.append({
                'Column': col,
                'Skewness': round(skew, 4),
                'Skew_Interpretation': skew_interpretation,
                'Kurtosis': round(kurt, 4),
                'Kurt_Interpretation': kurt_interpretation
            })
            
            print(f"{col:20s}: Skewness = {skew:7.4f} [{skew_interpretation:20s}], "
                  f"Kurtosis = {kurt:7.4f} [{kurt_interpretation}]")
        
        return pd.DataFrame(results)
    
    
    def multicollinearity_detection(self, threshold=0.9):
        """
        Detect multicollinearity using correlation threshold
        
        Parameters:
        -----------
        threshold : float
            Correlation threshold for multicollinearity (default: 0.9)
        """
        print("\n" + "="*80)
        print(f"MULTICOLLINEARITY DETECTION (threshold = {threshold})")
        print("="*80)
        
        corr_matrix = self.df[self.numeric_columns].corr().abs()
        high_corr_pairs = []
        
        for i in range(len(corr_matrix.columns)):
            for j in range(i+1, len(corr_matrix.columns)):
                if corr_matrix.iloc[i, j] >= threshold:
                    high_corr_pairs.append({
                        'Variable_1': corr_matrix.columns[i],
                        'Variable_2': corr_matrix.columns[j],
                        'Correlation': corr_matrix.iloc[i, j]
                    })
        
        if high_corr_pairs:
            print("\n⚠ High multicollinearity detected:")
            for pair in high_corr_pairs:
                print(f"  {pair['Variable_1']} <-> {pair['Variable_2']}: r = {pair['Correlation']:.4f}")
        else:
            print("✓ No significant multicollinearity detected")
        
        return pd.DataFrame(high_corr_pairs)
    
    
    def advanced_normality_tests(self):
        """
        Perform multiple normality tests (Shapiro-Wilk, Jarque-Bera, D'Agostino)
        """
        print("\n" + "="*80)
        print("ADVANCED NORMALITY TESTS")
        print("="*80)
        print("Tests: Shapiro-Wilk, Jarque-Bera, D'Agostino's Normality Test")
        print("-"*80)
        
        results = []
        for col in self.numeric_columns:
            data = self.df[col].dropna()
            
            # Shapiro-Wilk (for n < 5000)
            if len(data) < 5000:
                sw_stat, sw_p = stats.shapiro(data)
            else:
                sw_stat, sw_p = np.nan, np.nan
            
            # Jarque-Bera
            jb_stat, jb_p = jarque_bera(data)
            
            # D'Agostino
            da_stat, da_p = normaltest(data)
            
            # Consensus
            p_values = [p for p in [sw_p, jb_p, da_p] if not np.isnan(p)]
            consensus = "Normal" if all(p > 0.05 for p in p_values) else "Not Normal"
            
            results.append({
                'Column': col,
                'Shapiro_Wilk_p': round(sw_p, 4) if not np.isnan(sw_p) else 'N/A',
                'Jarque_Bera_p': round(jb_p, 4),
                'D_Agostino_p': round(da_p, 4),
                'Consensus': consensus
            })
            
            print(f"\n{col}:")
            if not np.isnan(sw_p):
                print(f"  Shapiro-Wilk: p = {sw_p:.4f}")
            print(f"  Jarque-Bera: p = {jb_p:.4f}")
            print(f"  D'Agostino: p = {da_p:.4f}")
            print(f"  → {consensus}")
        
        return pd.DataFrame(results)
    
    
    def plot_qq_plots(self):
        """
        Generate Q-Q plots for normality assessment
        """
        n = len(self.numeric_columns)
        rows = (n + 2) // 3
        plt.figure(figsize=(18, 5*rows))
        
        for i, col in enumerate(self.numeric_columns, 1):
            plt.subplot(rows, 3, i)
            stats.probplot(self.df[col].dropna(), dist="norm", plot=plt)
            plt.title(f"Q-Q Plot: {col}")
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig("qq_plots.png", dpi=300)
        print("✓ Saved: qq_plots.png")
        plt.show()
    
    
    def correlation_matrix_advanced(self):
        """
        Calculate Pearson, Spearman, and Kendall correlations
        """
        print("\n" + "="*80)
        print("ADVANCED CORRELATION ANALYSIS")
        print("="*80)
        
        results = {}
        
        for i, col1 in enumerate(self.numeric_columns):
            for col2 in self.numeric_columns[i+1:]:
                data1 = self.df[col1].dropna()
                data2 = self.df[col2].dropna()
                
                # Align indices
                common_idx = data1.index.intersection(data2.index)
                if len(common_idx) < 3:
                    continue
                
                d1 = data1.loc[common_idx]
                d2 = data2.loc[common_idx]
                
                pearson_r, pearson_p = pearsonr(d1, d2)
                spearman_r, spearman_p = spearmanr(d1, d2)
                kendall_r, kendall_p = kendalltau(d1, d2)
                
                results[f"{col1} vs {col2}"] = {
                    'Pearson_r': round(pearson_r, 4),
                    'Pearson_p': round(pearson_p, 4),
                    'Spearman_r': round(spearman_r, 4),
                    'Spearman_p': round(spearman_p, 4),
                    'Kendall_r': round(kendall_r, 4),
                    'Kendall_p': round(kendall_p, 4)
                }
        
        corr_df = pd.DataFrame(results).T
        print(corr_df.to_string())
        return corr_df
    
    
    def scale_features(self, method='standard', columns=None):
        """
        Scale/normalize features
        
        Parameters:
        -----------
        method : str
            'standard' (Z-score), 'minmax', 'robust'
        columns : list, optional
            Columns to scale. If None, scales all numeric columns
        """
        if columns is None:
            columns = self.numeric_columns
        
        print(f"\nScaling features using '{method}' method...")
        
        if method == 'standard':
            scaler = StandardScaler()
        elif method == 'minmax':
            scaler = MinMaxScaler()
        elif method == 'robust':
            scaler = RobustScaler()
        else:
            raise ValueError("Method must be 'standard', 'minmax', or 'robust'")
        
        df_scaled = self.df.copy()
        df_scaled[columns] = scaler.fit_transform(self.df[columns])
        
        self.df = df_scaled
        print(f"✓ Scaling complete for columns: {columns}")
        return self.df
    
    
    def data_quality_report(self):
        """
        Generate comprehensive data quality report
        """
        print("\n" + "="*80)
        print("DATA QUALITY REPORT")
        print("="*80)
        
        report = {}
        
        for col in self.numeric_columns:
            data = self.df[col].dropna()
            
            report[col] = {
                'Total_Values': len(self.df[col]),
                'Missing_Values': self.df[col].isnull().sum(),
                'Missing_Percentage': round((self.df[col].isnull().sum() / len(self.df[col])) * 100, 2),
                'Unique_Values': self.df[col].nunique(),
                'Zero_Values': (self.df[col] == 0).sum(),
                'Negative_Values': (self.df[col] < 0).sum(),
                'Mean': round(data.mean(), 4) if len(data) > 0 else np.nan,
                'Std': round(data.std(), 4) if len(data) > 0 else np.nan,
                'CV': round((data.std() / data.mean()) * 100, 2) if len(data) > 0 and data.mean() != 0 else np.nan,
                'Skewness': round(stats.skew(data), 4) if len(data) > 0 else np.nan,
                'Kurtosis': round(stats.kurtosis(data), 4) if len(data) > 0 else np.nan
            }
        
        report_df = pd.DataFrame(report).T
        print(report_df.to_string())
        
        # Save to CSV
        report_df.to_csv("data_quality_report.csv")
        print("\n✓ Saved: data_quality_report.csv")
        
        return report_df
    
    
    def plot_violin_plots(self):
        """
        Generate violin plots for distribution comparison
        """
        n = len(self.numeric_columns)
        rows = (n + 2) // 3
        plt.figure(figsize=(18, 5*rows))
        
        for i, col in enumerate(self.numeric_columns, 1):
            plt.subplot(rows, 3, i)
            sns.violinplot(y=self.df[col], inner="box")
            plt.title(f"{col} Distribution (Violin Plot)")
        
        plt.tight_layout()
        plt.savefig("violin_plots.png", dpi=300)
        print("✓ Saved: violin_plots.png")
        plt.show()
    
    
    def plot_pairplot(self, sample_size=1000):
        """
        Generate pairplot for all numeric variables
        
        Parameters:
        -----------
        sample_size : int
            Number of samples to plot (for large datasets)
        """
        data = self.df[self.numeric_columns]
        
        if len(data) > sample_size:
            data = data.sample(n=sample_size, random_state=42)
            print(f"\nSampling {sample_size} rows for pairplot...")
        
        sns.pairplot(data, diag_kind='kde', plot_kws={'alpha': 0.6})
        plt.suptitle('Pairwise Relationships', y=1.02, fontsize=16)
        plt.tight_layout()
        plt.savefig("pairplot.png", dpi=300)
        print("✓ Saved: pairplot.png")
        plt.show()
    
    
    def export_summary_report(self, filename=None):
        """
        Export comprehensive summary report to Excel
        
        Parameters:
        -----------
        filename : str, optional
            Output filename. If None, generates timestamped filename
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"eda_summary_report_{timestamp}.xlsx"
        
        print(f"\nGenerating summary report: {filename}")
        
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            # Basic statistics
            self.df[self.numeric_columns].describe().to_excel(writer, sheet_name='Basic_Statistics')
            
            # Correlation matrix
            self.df[self.numeric_columns].corr().to_excel(writer, sheet_name='Correlation_Matrix')
            
            # Missing data
            missing_df = pd.DataFrame({
                'Column': self.numeric_columns,
                'Missing_Count': [self.df[col].isnull().sum() for col in self.numeric_columns],
                'Missing_Percentage': [(self.df[col].isnull().sum() / len(self.df)) * 100 for col in self.numeric_columns]
            })
            missing_df.to_excel(writer, sheet_name='Missing_Data', index=False)
            
            # Skewness and Kurtosis
            skew_kurt = self.skewness_kurtosis_analysis()
            skew_kurt.to_excel(writer, sheet_name='Skewness_Kurtosis', index=False)
            
            # Data quality
            quality_report = self.data_quality_report()
            quality_report.to_excel(writer, sheet_name='Data_Quality')
        
        print(f"✓ Summary report saved: {filename}")
        return filename
    
    
    def statistical_summary_table(self):
        """
        Create comprehensive statistical summary table
        """
        print("\n" + "="*80)
        print("COMPREHENSIVE STATISTICAL SUMMARY")
        print("="*80)
        
        summary_data = []
        
        for col in self.numeric_columns:
            data = self.df[col].dropna()
            
            summary_data.append({
                'Column': col,
                'Count': len(data),
                'Mean': round(data.mean(), 4),
                'Median': round(data.median(), 4),
                'Std': round(data.std(), 4),
                'Min': round(data.min(), 4),
                'Max': round(data.max(), 4),
                'Q25': round(data.quantile(0.25), 4),
                'Q75': round(data.quantile(0.75), 4),
                'IQR': round(data.quantile(0.75) - data.quantile(0.25), 4),
                'Skewness': round(stats.skew(data), 4),
                'Kurtosis': round(stats.kurtosis(data), 4),
                'CV_%': round((data.std() / data.mean()) * 100, 2) if data.mean() != 0 else np.nan
            })
        
        summary_df = pd.DataFrame(summary_data)
        print(summary_df.to_string(index=False))
        return summary_df
    
    
    def run_full_analysis(self, target=None):
        """
        Run complete exploratory data analysis
        
        Parameters:
        -----------
        target : str, optional
            Target column for pairwise relationships. If None, uses first numeric column
        """
        if target is None:
            target = self.numeric_columns[0]
        
        print("\n" + "="*80)
        print("RUNNING FULL EXPLORATORY DATA ANALYSIS")
        print("="*80)
        self.basic_statistics()
        self.correlation_analysis()
        self.plot_distributions()
        self.boxplots()
        self.plot_correlation_heatmap()
        self.normality_tests()
        self.outlier_detection()
        self.plot_pairwise_relationships(target)