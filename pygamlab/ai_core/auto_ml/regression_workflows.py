import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.datasets import load_breast_cancer
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier


from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.inspection import permutation_importance
from sklearn.pipeline import Pipeline

from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, auc



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    confusion_matrix,
    precision_recall_curve,
    average_precision_score,
    roc_curve,
    auc,
    classification_report
)
from sklearn.inspection import permutation_importance
from sklearn.decomposition import PCA
#from sklearn.datasets import load_breast_cancer
from sklearn.datasets import fetch_california_housing
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score , mean_absolute_percentage_error
from sklearn.linear_model import LinearRegression, Ridge, Lasso, LogisticRegression
import joblib 
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression, Ridge, Lasso, LogisticRegression




class BASIC_regression:
    """
    General-purpose classification pipeline for Logistic Regression and other models.

    Workflow:
        1️ Initialize with X, y
        2️ split(test_size=0.2)
        3️ fitting(model='logistic')
        4️ predicting()
        5️ Accuracy(), plot_confusion_matrix(), plot_roc_curve_cv(), etc.
    """

    def __init__(self, X, y):
        if X is None or y is None:
            raise ValueError("Insert both X and y as input")

        self.X = X
        self.y = np.array(y).ravel()
        self.model = None
        self.ypred = None
        self.scaler = None
        self.df = None

    def split(self, test_size=0.2):
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            self.X, self.y, test_size=test_size, random_state=42
        )
        #self.scaler = MinMaxScaler()
        #self.X_train = self.scaler.fit_transform(self.X_train)
        #self.X_test = self.scaler.transform(self.X_test)
        print(f"Data split complete: train={len(self.X_train)}, test={len(self.X_test)}")


    def scaling(self,def_scaler):

        if def_scaler=='minmax':
            self.scaler = MinMaxScaler()
            self.X_train = self.scaler.fit_transform(self.X_train)
            self.X_test = self.scaler.transform(self.X_test)
            
        elif def_scaler==None:
            self.X_train =self.X_train 
            self.X_test=self.X_test
            

        else:
            self.scaler =def_scaler
            self.X_train = self.scaler.fit_transform(self.X_train)
            self.X_test = self.scaler.transform(self.X_test)



    def data_clear_check(self):
        df = pd.concat([pd.DataFrame(self.X), pd.Series(self.y, name="target")], axis=1)
        print("Null values:\n", df.isnull().sum())
        print("\nData Description:\n", df.describe())
        self.df = df

    def plot_correlation(self, save=False):
        if self.df is None:
            self.data_clear_check()
        corr = self.df.corr(numeric_only=True)
        plt.figure(figsize=(10, 10))
        sns.heatmap(np.abs(corr), cmap="RdYlGn", annot=True)
        plt.title("Absolute Pearson Correlation")
        if save:
            plt.savefig("correlation.jpg", format="jpg")
        plt.show()

    def plot_pair_plot(self, save=False):
        if self.df is None:
            self.data_clear_check()
        pairplot = sns.pairplot(self.df, hue="target", palette="coolwarm")
        if save:
            pairplot.savefig("pairplot.jpg", format="jpg")
        plt.show()

    

    def fitting(self, model):
        """
        Fit the model based on the string input.

        Parameters:
        - model: str, one of ['linear', 'ridge', 'lasso', 'logistic']
        - kwargs: additional arguments to pass to the model
        """
        self.model=model
        # Fit the model
        self.model.fit(self.X_train, self.y_train)

        # Make predictions
        self.ypred = self.model.predict(self.X_test)

        model_name = type(self.model).__name__
        print(f"{model_name} model fitted successfully.")


    def predicting(self,x_new):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
        prediction=self.model.predict(x_new)
        return prediction
        


    
    def get_score(self,metric='all'):
        if self.model is None:
            raise ValueError("Model not trained yet.")

        # Make sure test data exist
        if not hasattr(self, 'X_test') or not hasattr(self, 'y_test'):
            raise ValueError("X_test and y_test are missing. Please split your data before evaluation.")
        
        # Predictions
        y_pred = self.model.predict(self.X_test)

        # Metrics
        mae = mean_absolute_error(self.y_test, y_pred)

        if metric == 'mae':
            return mae
        
        mse = mean_squared_error(self.y_test, y_pred)

        if metric == 'mse':
            return mse
        
        rmse = np.sqrt(mse)

        if metric == 'rmse':
            return rmse
        
        
        mape=mean_absolute_percentage_error(self.y_test,y_pred)
        
        if metric=='mape':
            return mape
        

        r2 = r2_score(self.y_test, y_pred)

        if metric == 'r2':
            return r2

        print(f"Model Evaluation Metrics:")
        print(f"  Mean Absolute Error (MAE):  {mae:.4f}")
        print(f"  Mean Squared Error (MSE):   {mse:.4f}")
        print(f"  Root Mean Squared Error (RMSE): {rmse:.4f}")
        print(f"Mean Absolute Percentage Error(MAPE) : {mape:.4f}")
        print(f"  R² Score: {r2:.4f}")

        # Optionally return as a dictionary
        return {
            "MAE": mae,
            "MSE": mse,
            "RMSE": rmse,
            'MAPE':mape,
            "R2": r2
        }


    

    def intercept(self):
        if not hasattr(self.model, "coef_"):
            raise AttributeError("Model does not have coefficients.")
        print("Intercept:", self.model.intercept_)
        print("Coefficients:", self.model.coef_)
        
    

    def plot_true_vs_predict(self, name_model="Model", save=False):
        if self.ypred is None:
            raise ValueError("Predictions not found. Call predicting() first.")
        
        plt.figure(figsize=(8, 6))

        # Scatter plot: True vs Predicted
        plt.scatter(self.y_test, self.ypred, alpha=0.7, color='blue')
        plt.plot([self.y_test.min(), self.y_test.max()],
                [self.y_test.min(), self.y_test.max()], 'r--', lw=2)  # diagonal line
        plt.xlabel("True Values")
        plt.ylabel("Predicted Values")
        plt.title(f"True vs Predicted ({name_model})")
        plt.grid(True)

        if save:
            plt.savefig("true_vs_predict.jpg", format="jpg")
        plt.show()
        
        
        
        





class regression_CV:

    def __init__(self, X, y):
        if X is None or y is None:
            raise ValueError("Insert both X and y as input")

        self.X = X
        self.y = np.array(y).ravel()

        self.model = None
        self.ypred = None
        self.best_estimator = None
        self.best_score = None
        self.gs = None
        self.myparams = None
        self.df = None

    def data_clear_check(self):
        self.df = pd.concat([pd.DataFrame(self.X), pd.Series(self.y, name='target')], axis=1)

        print('-------------------')
        print(self.df.isnull().sum())
        print('-------------------')
        print(self.df.isna().sum())
        print('-------------------')
        print(self.df.describe())

    def plot_correlation(self, save=False):
        if self.df is None:
            self.data_clear_check()
        yX_corr = self.df.corr(method='pearson')
        yX_abs_corr = np.abs(yX_corr)

        plt.figure(figsize=(10, 10))
        sns.heatmap(yX_abs_corr, cmap='RdYlGn', square=True)
        plt.title('Pearson Correlation Heat Map (absolute values)', fontsize=15, fontweight='bold')
        if save:
            plt.savefig("correlation.jpg", format='jpg')
        plt.show()

    def plot_pair_plot(self, hue_col='target', save=False):
        if self.df is None:
            self.data_clear_check()
        pairplot = sns.pairplot(self.df, hue=hue_col, palette='coolwarm')
        if save:
            pairplot.savefig("pairplot.jpg", format='jpg')
        plt.show()

    def scaling(self, def_scaler='minmax'):
        if def_scaler == 'minmax':
            scaler = MinMaxScaler()
        else:
            scaler = def_scaler

        self.X = scaler.fit_transform(self.X)
        self.y = np.array(self.y).ravel()

    def fitting(self, model, myparams=None, cv=5, scoring='neg_mean_absolute_percentage_error'):
        self.model = model
        self.myparams = myparams if myparams is not None else {}

        self.gs = GridSearchCV(
            self.model,
            self.myparams,
            cv=cv,
            scoring=scoring,
            return_train_score=True,
            n_jobs=-1,
        )
        self.gs.fit(self.X, self.y)

        self.best_estimator = self.gs.best_estimator_
        self.best_score = self.gs.best_score_
        self.ypred = cross_val_predict(self.best_estimator, self.X, self.y, cv=cv)

        print("Regression model trained with GridSearchCV.")
        print(f"Best score: {self.best_score:.4f}")

    def best(self):
        if self.gs is None:
            raise ValueError("Model not trained yet. Call fitting() first.")

        self.best_estimator = self.gs.best_estimator_
        self.best_score = self.gs.best_score_
        return self.best_estimator, self.best_score

    def predicting(self, x_new):
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")
        prediction = self.best_estimator.predict(x_new)
        return prediction

    def get_score(self, metric='all'):
        if self.best_estimator is None or self.ypred is None:
            raise ValueError("Model not trained yet. Call fitting() first.")

        mae = mean_absolute_error(self.y, self.ypred)
        if metric == 'mae':
            return mae

        mse = mean_squared_error(self.y, self.ypred)
        if metric == 'mse':
            return mse

        rmse = np.sqrt(mse)
        if metric == 'rmse':
            return rmse

        mape = mean_absolute_percentage_error(self.y, self.ypred)
        if metric == 'mape':
            return mape

        r2 = r2_score(self.y, self.ypred)
        if metric == 'r2':
            return r2

        print("Model Evaluation Metrics:")
        print(f"  Mean Absolute Error (MAE):  {mae:.4f}")
        print(f"  Mean Squared Error (MSE):   {mse:.4f}")
        print(f"  Root Mean Squared Error (RMSE): {rmse:.4f}")
        print(f"  Mean Absolute Percentage Error (MAPE): {mape:.4f}")
        print(f"  R² Score: {r2:.4f}")

        return {
            "MAE": mae,
            "MSE": mse,
            "RMSE": rmse,
            "MAPE": mape,
            "R2": r2,
        }

    def plot_grid_search(self, save=False):
        if self.gs is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
        if self.myparams is None or len(self.myparams.keys()) != 2:
            raise ValueError("Heatmap only supported when tuning exactly 2 parameters.")

        key_params = list(self.myparams.keys())
        grid_params_x = self.myparams[key_params[0]]
        grid_params_y = self.myparams[key_params[1]]

        df_results = pd.DataFrame(self.gs.cv_results_)
        ar_scores = np.array(df_results['mean_test_score']).reshape(len(grid_params_y), len(grid_params_x))

        plt.figure(figsize=(8, 6))
        sns.heatmap(ar_scores, annot=True, fmt='.3f', cmap='coolwarm',
                    xticklabels=grid_params_x, yticklabels=grid_params_y)
        plt.title('Grid Search Heatmap (Regression)')
        plt.xlabel(key_params[0])
        plt.ylabel(key_params[1])
        plt.tight_layout()

        if save:
            plt.savefig("grid_search_heatmap_regression.jpg", format='jpg')
        plt.show()

    def intercept(self):
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")
        if not hasattr(self.best_estimator, "coef_") or not hasattr(self.best_estimator, "intercept_"):
            raise AttributeError("This regression model does not have coef_ or intercept_ attributes.")

        print("Intercept:")
        print(self.best_estimator.intercept_)
        print("\nCoefficients:")
        print(self.best_estimator.coef_)

    def plot_true_vs_predict(self, name_model="Model", save=False):
        if self.ypred is None:
            raise ValueError("Predictions not found. Call fitting() first.")

        plt.figure(figsize=(8, 6))
        plt.scatter(self.y, self.ypred, alpha=0.7, color='blue')
        plt.plot([self.y.min(), self.y.max()], [self.y.min(), self.y.max()], 'r--', lw=2)
        plt.xlabel("True Values")
        plt.ylabel("Predicted Values")
        plt.title(f"True vs Predicted ({name_model})")
        plt.grid(True)

        if save:
            plt.savefig("true_vs_predict.jpg", format="jpg")
        plt.show()

    def plot_permutation_importance(self, save=False):
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")

        result = permutation_importance(
            self.best_estimator, self.X, self.y, n_repeats=10, random_state=42, n_jobs=2
        )
        perm_sorted_idx = result.importances_mean.argsort()

        plt.figure(figsize=(8, 6))
        plt.boxplot(result.importances[perm_sorted_idx].T, vert=False)
        plt.yticks(range(len(perm_sorted_idx)), pd.DataFrame(self.X).columns[perm_sorted_idx])
        plt.axvline(x=0, color="k", linestyle="--")
        plt.xlabel("Permutation Importance (Decrease in Score)")
        plt.title("Permutation Feature Importance (Regression)")
        plt.tight_layout()

        if save:
            plt.savefig("permutation_importance_regression.jpg", format='jpg')
        plt.show()


class regression_PIPE:

    def __init__(self, X, y):
        if X is None or y is None:
            raise ValueError("Insert both X and y as input")

        self.X = X
        self.y = np.array(y).ravel()

        self.model = None
        self.ypred = None
        self.best_estimator = None
        self.best_score = None
        self.gs = None
        self.myparams = None
        self.pipe = None
        self.cv_folds = None


    def data_clear_check(self):
        self.df = pd.concat([pd.DataFrame(self.X), pd.Series(self.y, name='target')], axis=1)

        print('-------------------')
        print(self.df.isnull().sum())
        print('-------------------')
        print(self.df.isna().sum())
        print('-------------------')
        print(self.df.describe())


    def plot_correlation(self, save=False):
        yX_corr = self.df.corr(method='pearson')
        yX_abs_corr = np.abs(yX_corr)

        plt.figure(figsize=(10, 10))
        sns.heatmap(yX_abs_corr, cmap='RdYlGn', square=True)
        plt.title('Pearson Correlation Heat Map (absolute values)', fontsize=15, fontweight='bold')
        if save:
            plt.savefig("correlation.jpg", format='jpg')
        plt.show()


    def plot_pair_plot(self, hue_col='target', save=False):
        pairplot = sns.pairplot(self.df, hue=hue_col, palette='coolwarm')
        if save:
            pairplot.savefig("pairplot.jpg", format='jpg')
        plt.show()


    def make_pipe(self, model,poly=None,specific_preprocessor=None):

        if poly!=None and specific_preprocessor==None:
            self.model = model
            self.poly=PolynomialFeatures()
            self.pipe = Pipeline([
                ("scaler", MinMaxScaler()),
                ('poly', self.poly),
                ("model", self.model)
            ])


        elif poly==None and specific_preprocessor!=None:
            self.model = model
            self.scaler=specific_preprocessor
            self.pipe = Pipeline([
                ("scaler", self.scaler),
                ("model", self.model)
            ])


        elif poly!=None and specific_preprocessor!=None:
            self.model = model
            self.scaler=specific_preprocessor
            self.poly=PolynomialFeatures()

            self.pipe = Pipeline([
                ("scaler", self.scaler),
                ("poly", self.poly),
                ("model", self.model)
            ])

        else:
            """Create pipeline with MinMaxScaler and model."""
            self.model = model
            self.pipe = Pipeline([
                ("scaler", MinMaxScaler()),
                ("model", self.model)
            ])

        print("Pipeline is constructed successfully.")


    
    def fitting(self, myparams=None, cv=5,scoring='neg_mean_absolute_percentage_error'):
        """
        Trains a regression model with GridSearchCV.
        Uses R² score by default as the evaluation metric.
        """
        # ✅ Default model
        
        # ✅ Default hyperparameters (empty if none specified)
        self.myparams = myparams if myparams is not None else {}

        # ✅ Store CV folds for later use in evaluation
        self.cv_folds = cv

        # ✅ Run GridSearchCV
        self.gs = GridSearchCV(
            estimator=self.pipe,
            param_grid=self.myparams,
            cv=cv,
            scoring=scoring,
            return_train_score=True,
            n_jobs=-1
        )

        self.gs.fit(self.X, self.y)

        # ✅ Save best model and score
        self.best_estimator = self.gs.best_estimator_
        self.best_score = self.gs.best_score_

        self.ypred = cross_val_predict(self.best_estimator, self.X, self.y)

        print("✅ Regression model trained with GridSearchCV.")
        print(f"Best MAPE Score: {self.best_score:.4f}")



    def best(self):
        """Return best estimator and score."""
        if self.gs is None:
            raise ValueError("Model not trained yet. Call fitting() first.")

        self.best_estimator = self.gs.best_estimator_
        self.best_score = self.gs.best_score_
        return self.best_estimator, self.best_score


    def predicting(self,x_new):
        """Cross-validated prediction using best model."""
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")
        prediction=self.best_estimator.predict(x_new)
        return prediction

    def evaluate_cv_train_test_metrics(self, shuffle=False, random_state=42):
        """
        Calculate MAPE and MAE for both training and test sets using cross-validation.
        Uses the same CV strategy as GridSearchCV to get metrics from CV results.
        
        Parameters:
        -----------
        shuffle : bool, default=False
            Whether to shuffle data before splitting (matches GridSearchCV default)
        random_state : int, default=42
            Random state for reproducibility of CV splits
        
        Returns:
        --------
        dict : Dictionary containing MAPE and MAE for train and test sets from CV
        """
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() or fitting() first.")
        
        if self.cv_folds is None:
            raise ValueError("CV folds not found. Call fitting() first.")
        
        # Create KFold object with same number of folds as used in GridSearchCV
        # Default shuffle=False to match GridSearchCV behavior when cv is an integer
        kfold = KFold(n_splits=self.cv_folds, shuffle=shuffle, random_state=random_state if shuffle else None)
        
        # Lists to store metrics for each fold
        train_mape_list = []
        train_mae_list = []
        test_mape_list = []
        test_mae_list = []
        
        # Perform cross-validation manually to get train and test predictions
        for train_idx, test_idx in kfold.split(self.X):
            X_train_fold = self.X[train_idx]
            X_test_fold = self.X[test_idx]
            y_train_fold = self.y[train_idx]
            y_test_fold = self.y[test_idx]
            
            # Clone and fit the best estimator on training fold
            from sklearn.base import clone
            model_fold = clone(self.best_estimator)
            model_fold.fit(X_train_fold, y_train_fold)
            
            # Predict on both training and test sets
            y_train_pred = model_fold.predict(X_train_fold)
            y_test_pred = model_fold.predict(X_test_fold)
            
            # Calculate metrics for training fold
            train_mape = mean_absolute_percentage_error(y_train_fold, y_train_pred)
            train_mae = mean_absolute_error(y_train_fold, y_train_pred)
            
            # Calculate metrics for test fold
            test_mape = mean_absolute_percentage_error(y_test_fold, y_test_pred)
            test_mae = mean_absolute_error(y_test_fold, y_test_pred)
            
            # Append to lists
            train_mape_list.append(train_mape)
            train_mae_list.append(train_mae)
            test_mape_list.append(test_mape)
            test_mae_list.append(test_mae)
        
        # Calculate mean across all folds
        mean_train_mape = np.mean(train_mape_list)
        mean_train_mae = np.mean(train_mae_list)
        mean_test_mape = np.mean(test_mape_list)
        mean_test_mae = np.mean(test_mae_list)
        
        # Calculate std across all folds (optional, for showing variability)
        std_train_mape = np.std(train_mape_list)
        std_train_mae = np.std(train_mae_list)
        std_test_mape = np.std(test_mape_list)
        std_test_mae = np.std(test_mae_list)
        
        # Print results
        print("=" * 70)
        print("CROSS-VALIDATION METRICS (from CV results)")
        print("=" * 70)
        print("TRAINING SET METRICS (averaged across CV folds):")
        print(f"  MAPE (Mean Absolute Percentage Error): {mean_train_mape:.4f} ± {std_train_mape:.4f}")
        print(f"  MAE (Mean Absolute Error):              {mean_train_mae:.4f} ± {std_train_mae:.4f}")
        print("-" * 70)
        print("TEST SET METRICS (averaged across CV folds):")
        print(f"  MAPE (Mean Absolute Percentage Error): {mean_test_mape:.4f} ± {std_test_mape:.4f}")
        print(f"  MAE (Mean Absolute Error):              {mean_test_mae:.4f} ± {std_test_mae:.4f}")
        print("=" * 70)
        
        # Return as dictionary
        metrics = {
            'train': {
                'MAPE_mean': mean_train_mape,
                'MAPE_std': std_train_mape,
                'MAE_mean': mean_train_mae,
                'MAE_std': std_train_mae,
                'MAPE_all_folds': train_mape_list,
                'MAE_all_folds': train_mae_list
            },
            'test': {
                'MAPE_mean': mean_test_mape,
                'MAPE_std': std_test_mape,
                'MAE_mean': mean_test_mae,
                'MAE_std': std_test_mae,
                'MAPE_all_folds': test_mape_list,
                'MAE_all_folds': test_mae_list
            }
        }
        
        return metrics


    def plot_grid_search(self, save=False):
        """
        Plot a heatmap of GridSearchCV results for regression models.
        Works only if exactly 2 hyperparameters were tuned.
        """
        if self.gs is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
        if len(self.myparams.keys()) != 2:
            raise ValueError("Heatmap only supported when tuning exactly 2 parameters.")

        # Extract parameter names and grid values
        key_params = list(self.myparams.keys())
        grid_params_x = self.myparams[key_params[0]]
        grid_params_y = self.myparams[key_params[1]]

        # Extract GridSearch results into DataFrame
        df_results = pd.DataFrame(self.gs.cv_results_)

        # Reshape mean test scores to match parameter grid shape
        ar_scores = np.array(df_results['mean_test_score']).reshape(
            len(grid_params_y), len(grid_params_x)
        )

        # Plot heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(ar_scores, annot=True, fmt='.3f', cmap='coolwarm',
                    xticklabels=grid_params_x, yticklabels=grid_params_y)
        plt.title('Grid Search Heatmap (R² Score)')
        plt.xlabel(key_params[0])
        plt.ylabel(key_params[1])
        plt.tight_layout()

        if save:
            plt.savefig("grid_search_heatmap_regression.jpg", format='jpg')
        plt.show()



    


    def intercept(self):
        """Print regression model coefficients and intercept."""
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")

        model = self.best_estimator

        # Ensure model has coefficients
        if not hasattr(model, 'coef_') or not hasattr(model, 'intercept_'):
            raise AttributeError("This regression model does not have coef_ or intercept_ attributes.")

        print("Intercept:")
        print(model.intercept_)
        print("\nCoefficients:")
        print(model.coef_)



    def plot_true_vs_predict(self, name_model="Model", save=False ,professional_plotting=False):
        """
        Plot True vs Predicted values for regression (using full dataset).
        Works with cross_val_predict results (no explicit y_test needed).
        """
        if self.ypred is None:
            raise ValueError("Predictions not found. Call predicting() first.")
        if self.y is None:
            raise ValueError("True values not found. Check initialization or fitting.")
            
        if professional_plotting ==True:
            y_true = np.array(self.y).reshape(-1)
            y_pred = np.array(self.ypred).reshape(-1)

            # Metrics
            r2 = r2_score(y_true, y_pred)
            mae = mean_absolute_error(y_true, y_pred)
            mape= mean_absolute_percentage_error(y_true, y_pred)

            # Figure
            plt.figure(figsize=(9, 7))
            plt.style.use("seaborn-v0_8")

            # Scatter plot
            plt.scatter(y_true, y_pred, alpha=0.65, edgecolor="k", linewidth=0.5)

            # Perfect prediction line
            min_val = min(y_true.min(), y_pred.min())
            max_val = max(y_true.max(), y_pred.max())
            plt.plot([min_val, max_val], [min_val, max_val], "r--", lw=2, label="Ideal Prediction")

            # Regression line
            m, b = np.polyfit(y_true, y_pred, 1)
            plt.plot(y_true, m * y_true + b, color="green", lw=2, label="Trend Line")

            # Labels & Title
            plt.xlabel("True Values", fontsize=12)
            plt.ylabel("Predicted Values", fontsize=12)
            plt.title(f"True vs Predicted — {name_model}", fontsize=14, fontweight="bold")

            # Add metrics text box
            text_box = (
                f"MAE: {mae:.3f}\n"
                f"MAPE: {mape:.3f}"
            )
            plt.gca().text(
                0.05, 0.95, text_box, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
            )

            plt.legend()
            plt.grid(True, alpha=0.3)

            # Save plot if requested
            if save:
                file_name = f"true_vs_predicted_{name_model.replace(' ', '_')}.jpg"
                plt.savefig(file_name, format="jpg", dpi=300, bbox_inches="tight")

            plt.show()
            
            
            
        else:

            plt.figure(figsize=(8, 6))
    
            # Scatter plot: True vs Predicted
            plt.scatter(self.y, self.ypred, alpha=0.7, color='blue')
            plt.plot([self.y.min(), self.y.max()],
                    [self.y.min(), self.y.max()], 'r--', lw=2)  # diagonal line
            plt.xlabel("True Values")
            plt.ylabel("Predicted Values")
            plt.title(f"True vs Predicted ({name_model})")
            plt.grid(True)
    
            if save:
                plt.savefig("true_vs_predict.jpg", format="jpg")
            plt.show()
        
   
    

    def plot_permutation_importance(self, save=False):
        """Plot permutation feature importance for regression models."""
        if self.best_estimator is None:
            raise ValueError("Best estimator not found. Call best() first.")

        model = self.best_estimator
        result = permutation_importance(
            model, self.X, self.y, n_repeats=10, random_state=42, n_jobs=2
        )
        perm_sorted_idx = result.importances_mean.argsort()

        plt.figure(figsize=(8, 6))
        plt.boxplot(result.importances[perm_sorted_idx].T, vert=False)
        plt.yticks(range(len(perm_sorted_idx)), pd.DataFrame(self.X).columns[perm_sorted_idx])
        plt.axvline(x=0, color="k", linestyle="--")
        plt.xlabel("Permutation Importance (Decrease in R² Score)")
        plt.title("Permutation Feature Importance (Regression)")
        plt.tight_layout()
        
        if save:
            plt.savefig("permutation_importance_regression.jpg", format='jpg')
            
        plt.show()
        
    #it needs to be change and df 
        
        
    def plot_3d_prediction_maps(self, df , feature_x, feature_y, feature_target, grid_size=40, save=False):
        """
        Create 3D surface, 3D scatter and heatmap for two features.
        Works even if X and y were numpy arrays originally because df stores columns.
    
        Parameters
        ----------
        feature_x : str
            Name of feature on X-axis (must be column in df)
        feature_y : str
            Name of feature on Y-axis
        grid_size : int
            Number of interpolation points
        save : bool
            Save figure as JPG
        """
    
        if self.best_estimator is None:
            raise ValueError("Train model first using fitting().")
    
        
    
        # Extract selected features
        x_vals = df[feature_x].values
        y_vals = df[feature_y].values
        z_vals = df[feature_target].values
    
        # Build mesh grid
        x_lin = np.linspace(x_vals.min(), x_vals.max(), grid_size)
        y_lin = np.linspace(y_vals.min(), y_vals.max(), grid_size)
    
        xx, yy = np.meshgrid(x_lin, y_lin)
    
        # Prepare full X input for prediction
        X_grid = pd.DataFrame(
            np.zeros((grid_size * grid_size, df.shape[1] - 1)),
            columns=df.drop(columns=[feature_target]).columns
        )
    
        # Fill selected two columns with mesh
        X_grid[feature_x] = xx.flatten()
        X_grid[feature_y] = yy.flatten()
    
        # For other columns, fill with their mean
        for col in X_grid.columns:
            if col not in [feature_x, feature_y]:
                X_grid[col] = df[col].mean()
    
        # Predict on entire grid
        z_pred = self.best_estimator.predict(X_grid)
        zz = z_pred.reshape(grid_size, grid_size)
    
        # --- PLOTTING ---
        fig = plt.figure(figsize=(20, 6))
    
        # 1) 3D prediction surface
        #ax1 = fig.add_subplot(1, 3, 1, projection="3d")
        ax1 = fig.add_subplot(1, 2, 1, projection="3d")
        surf = ax1.plot_surface(xx, yy, zz, cmap="plasma", alpha=0.85)
        ax1.set_xlabel(feature_x)
        ax1.set_ylabel(feature_y)
        ax1.set_zlabel("Predicted Y")
        ax1.set_title("3D Prediction Surface", fontsize=14)
        fig.colorbar(surf, ax=ax1, shrink=0.6)
    
        # 2) 3D scatter of real data
        '''
        ax2 = fig.add_subplot(1, 3, 2, projection="3d")
        sc = ax2.scatter(x_vals, y_vals, z_vals, c=z_vals, cmap="viridis", s=40, edgecolor="k")
        ax2.set_xlabel(feature_x)
        ax2.set_ylabel(feature_y)
        ax2.set_zlabel("True Y")
        ax2.set_title("3D Scatter (Real Data)", fontsize=14)
        fig.colorbar(sc, ax=ax2, shrink=0.6)
        '''
    
        # 3) 2D heatmap
        #ax3 = fig.add_subplot(1, 3, 3)
        ax3 = fig.add_subplot(1, 2, 2)
        hm = ax3.contourf(xx, yy, zz, levels=40, cmap="plasma")
        ax3.set_xlabel(feature_x)
        ax3.set_ylabel(feature_y)
        ax3.set_title("Prediction Heatmap", fontsize=14)
        fig.colorbar(hm, ax=ax3)
    
        plt.tight_layout()
    
        if save:
            plt.savefig("3d_prediction_maps.jpg", dpi=300, bbox_inches="tight")
    
        plt.show()


