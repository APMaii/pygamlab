import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.datasets import load_breast_cancer
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PolynomialFeatures
from sklearn.decomposition import PCA

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier


from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.inspection import permutation_importance
from sklearn.pipeline import Pipeline

from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, auc



class BASIC_classification:
    
    '''
    - - Classification with Train test split:
        first use split() and then use the other attributes.
        Input--> X as iput, y as output
        
        Attributes:
            
            split():
                
            data_clear_check():
                
                
            fitting():
                
            predicting():
                
            accuracy():
                
            data_clear_check():
                
            plot_correlation():
                
            plot_pair_plot():
                
            plot_confusion_matrix():
                
            plot_roc_curve_cv():
    
    
    
    '''

    def __init__(self, X, y):
        
        if X is None or y is None:
            raise ValueError("Insert Both x and y as input")
            
        self.X=X

        self.y= np.array(y).ravel()
        
        self.model=None
        self.ypred=None
        #self.model.intercept_ = None  
        #self.model.coef_ =None  
        #self.model.feature_importances_ = None
        

        
    def split(self,test_size):
        self.X_train, self.X_test, self.y_train, self.y_test=train_test_split(self.X, self.y, test_size=test_size, random_state=42)
        
        
        
    def scaling(self,def_scaler='minmax'):
        
        if def_scaler=='minmax':
            scaler = MinMaxScaler()
            self.X_train= scaler.fit_transform(self.X_train)
            self.X_test=scaler.transform(self.X_test)
            
        else:
            scaler=def_scaler
            self.X_train= scaler.fit_transform(self.X_train)
            self.X_test=scaler.transform(self.X_test)
            
            
        


    def data_clear_check(self):
       self.df=pd.concat([pd.DataFrame(self.X),self.y],axis=1)
       
       #df=self.df
       print('-------------------')
       self.df.isnull().sum()
       print('-------------------')
       self.df.isna().sum()
       print('-------------------')
       self.df.describe()
       
       
    def plot_correlation(self,save=False):
        yX_corr = self.df.corr(method='pearson')

          # Convert to abolute values
        yX_abs_corr = np.abs(yX_corr) 
          

        plt.figure(figsize=(10, 10))
        plt.imshow(yX_abs_corr, cmap='RdYlGn', interpolation='none', aspect='auto')
        plt.colorbar()
        plt.xticks(range(len(yX_abs_corr)), yX_abs_corr.columns, rotation='vertical')
        plt.yticks(range(len(yX_abs_corr)), yX_abs_corr.columns);
        plt.suptitle('Pearson Correlation Heat Map (absolute values)', fontsize=15, fontweight='bold')
        if save is True:
            plt.savefig("correlation.jpg", format='jpg')
        plt.show()
        
            
    def plot_pair_plot(self,save=False):
        pairplot = sns.pairplot(self.df, hue='diognosis', palette='coolwarm')
        if save is True:
            # Save the pairplot as jpg
            pairplot.savefig("pairplot.jpg", format='jpg')

        # Show the pairplot
        plt.show()

  
    def fitting(self,model):
        self.model=model
        self.model.fit(self.X_train,self.y_train)
        
        print('The model is trained and is ready now!')
        
        self.ypred=self.model.predict(self.X_test)
        
        
    def Accuracy(self):
        
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        self.score=self.model.score(self.X_test,self.y_test)

        return self.score
    
    
    def predicting(self,new_x):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        predicted=self.model.predict(new_x)
        
        return predicted
    
    

    def plot_confusion_matrix(self,save=False):
        
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
            
        cm = confusion_matrix(self.y_test, self.ypred)
        confusion = sns.heatmap(cm, annot = True,fmt='.4g')
        if save is True:
            # Save the pairplot as jpg
            confusion.savefig("confusion_matrix.jpg", format='jpg')

        # Show the pairplot
        plt.show()

         
    '''
   
        
      y_pred = gs_svc.predict(X).astype(int)
      y_score = gs_svc.predict_proba(X)
    
    plot_roc_and_precision_recall(Y, y_score[:,1])  # provide the column for the scores belonging only to the positive class

    '''

    def plot_precision_recall_cv(self,save=False):
        #some chekc  for models predict proba AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
        
        
        y_scores = self.model.predict_proba(self.X_test)
        
        precision, recall, _ = precision_recall_curve(self.y_test, y_scores[:, 1])
        average_precision = average_precision_score(self.y_test, y_scores[:, 1])

        plt.figure()
        plt.step(recall, precision, where='post', color='b', alpha=0.5)
        plt.fill_between(recall, precision, step='post', alpha=0.5, color='b')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("precision_recall_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()


    def plot_roc_curve_cv(self,save=False):
        #also some chekc for some models AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            

            
        y_scores = self.model.predict_proba(self.X_test)
        
        fpr, tpr, _ = roc_curve(self.y_test, y_scores[:, 1])
        roc_auc = auc(fpr, tpr)

        plt.figure()
        lw = 2
        plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend(loc="lower right")
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("roc_curve_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()



    def intercept(self):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.model.intercept_ is None or self.model.coef_ is None :
            raise AttributeError('Your model can not provide coefficients and intercept, Make sure you use appropriate model!')
        
            
        print("intercept_:")
        print(self.model.intercept_ )
        print()
        print("coef_:")
        print(self.model.coef_)
        
        
    def plot_feature_importance(self,save=False):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
        
        if self.model.feature_importances_ is None:
            raise AttributeError('This function is used only for Random Forest classifier!')
            
            
        dataset=load_breast_cancer()
        n_feature = 30
        plt.barh(range(n_feature), self.model.feature_importances_, align='center')
        plt.yticks(np.arange(n_feature), dataset.feature_names)
        plt.xlabel('Feature Importance')
        plt.ylabel('Feature')
        if save is True:
            plt.savefig("feature_importance.jpg", format='jpg')
        plt.show()
        
        
    
    def plot_coeff_weitgh(self,save=False):
        
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.model.coefs_ is None:
            raise AttributeError('This function is used only for MLP classifier!')
            
        dataset=load_breast_cancer()
        
        plt.figure(figsize=(20,5))
        plt.imshow(self.model.coefs_[0],interpolation='None',cmap='GnBu')
        plt.yticks(range(30), dataset.feature_names)
        plt.xlabel('Colums in weight matrix')
        plt.ylabel('Input feature')
        plt.colorbar()
        if save is True:
            plt.savefig("coeff_weight.jpg", format='jpg')
        plt.show()
        


    def plot_true_vs_predict(self,name_model,save=False):
        
        #for test 
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
        # Perform PCA for dimensionality reduction
        pca = PCA(n_components=2)
        self.X_train_pca = pca.fit_transform(self.X_train)
        self.X_test_pca=pca.fit(self.X_test)
        
        #self.model.fit(self.X_train_pca,self.y_train)
        #self.y_pred_pca=self.model.predict(self.X_test_pca)
        
 
        # Visualize the true labels versus predicted labels for training set
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.scatter(self.X_test_pca[:, 0], self.X_test_pca[:, 1], c=np.array(self.y_test), cmap='viridis', label='True Labels')
        plt.title(f'True Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.scatter(self.X_test_pca[:, 0], self.X_test_pca[:, 1], c=self.ypred, cmap='viridis', label='Predicted Labels')
        plt.title(f'Predicted Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.tight_layout()
        if save is True:
            plt.savefig("true_vs_predict.jpg", format='jpg')
        plt.show()
        


        
    def plot_permutation_importance(self,save=False):
        if self.model is None:
            raise ValueError("Best is not fined. Call fitting() first.")
        
        
        # Calculate Permutation Feature Importance
        result = permutation_importance(self.model, self.X, self.y, n_repeats=10, random_state=42, n_jobs=2)
        perm_sorted_idx = result.importances_mean.argsort()
        
        fig, ax = plt.subplots(figsize=(8, 6))  # Create a new figure and axes
        ax.boxplot(
            result.importances[perm_sorted_idx].T,
            vert=False,
            labels=pd.DataFrame(self.X).columns[perm_sorted_idx],
            )
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Permutation Importance")
        ax.set_title("Permutation Feature Importance")
        if save is True:
            plt.savefig("permutation_importance.jpg", format='jpg')
        plt.show()













class classification_CV:

    def __init__(self, X, y):
        
        if X is None or y is None:
            raise ValueError("Insert Both x and y as input")
            
        #self.X=X
        #self.y=y

        
        self.X=X
        self.y=y
        
        self.model=None
        self.ypred=None
        self.best_estimator=None
        self.gs=None
        self.myparams=None
        
        #self.best_estimator.intercept_ = None  
        #self.best_estimator.coef_ =None  
        #self.best_estimator.feature_importances_ = None
        

        
        
    def data_clear_check(self):
       self.df=pd.concat([pd.DataFrame(self.X),self.y],axis=1)
       
       #df=self.df
       print('-------------------')
       self.df.isnull().sum()
       print('-------------------')
       self.df.isna().sum()
       print('-------------------')
       self.df.describe()
       
       
    def plot_correlation(self,save=False):
        yX_corr = self.df.corr(method='pearson')

          # Convert to abolute values
        yX_abs_corr = np.abs(yX_corr) 
          

        plt.figure(figsize=(10, 10))
        plt.imshow(yX_abs_corr, cmap='RdYlGn', interpolation='none', aspect='auto')
        plt.colorbar()
        plt.xticks(range(len(yX_abs_corr)), yX_abs_corr.columns, rotation='vertical')
        plt.yticks(range(len(yX_abs_corr)), yX_abs_corr.columns);
        plt.suptitle('Pearson Correlation Heat Map (absolute values)', fontsize=15, fontweight='bold')
        if save is True:
            plt.savefig("correlation.jpg", format='jpg')
        plt.show()
        
            
    def plot_pair_plot(self,save=False):
        pairplot = sns.pairplot(self.df, hue='diognosis', palette='coolwarm')
        if save is True:
            # Save the pairplot as jpg
            pairplot.savefig("pairplot.jpg", format='jpg')

        # Show the pairplot
        plt.show()
        
        
    def scaling(self,def_scaler='minmax'):
        
        if def_scaler=='minmax':
            scaler = MinMaxScaler()
            scaler.fit(X)
            self.X = scaler.transform(X)
            
            self.y= np.array(y).ravel()
            
        else:
            
            scaler=def_scaler
            scaler.fit(X)
            self.X = scaler.transform(X)
            
            self.y= np.array(y).ravel()
            
        

  
    def fitting(self,model,myparams):
        self.model=model
        self.myparams=myparams
        
        self.gs= GridSearchCV(self.model, self.myparams , return_train_score=True)
        self.gs.fit(self.X,self.y)
        
        print('The model is trained and is ready now!')
        
        self.best_estimator=self.gs.best_estimator_
        self.best_score=self.gs.best_score_
        
        self.ypred=cross_val_predict(self.best_estimator,self.X,self.y)
        
        
    def best(self):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        self.best_estimator=self.gs.best_estimator_
        self.best_score=self.gs.best_score_
        
        return self.best_estimator,self.best_score
    
    
    def predicting(self,x_new):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
               
            
        prediction=self.gs.best_estimator_.predict(x_new)
        return prediction

    
    def plot_grid_search(self, save=False):
      '''
      Plot 2d grid search heatmap
      
      Parameters:
      grid_search: instance of sklearn's GridSearchCV
      grid_params: dictionary of grid search parameters  
      x_param: name of x-axis parameter in grid_params
      y_param: name of y-axis parameter in grid_params
      is_verbose (optional): print results
      
      Return:
      grid_search.best_score_: best score found
      grid_search.best_estimator_: best estimator found
      '''

      key_params=list(self.myparams.keys())
      
      grid_params_x = self.myparams[key_params[0]]
      grid_params_y = self.myparams[key_params[1]]
      
      if self.gs is None:
          raise ValueError("Model not trained yet. Call fitting() first.")
          
      if self.myparams is None:
          raise ValueError("Model not trained yet. Call fitting() first.")
      
        
      df_results = pd.DataFrame(self.gs.cv_results_)
      ar_scores = np.array(df_results.mean_test_score).reshape(len(grid_params_y), len(grid_params_x))
      sns.heatmap(ar_scores, annot=True, fmt='.3f', xticklabels=grid_params_x, yticklabels=grid_params_y)
      print()
      plt.suptitle('Grid search heatmap')
      plt.xlabel(key_params[0])
      plt.ylabel(key_params[1])
      if save is True:
          # Save the pairplot as jpg
          plt.savefig("grid_search_heat_map.jpg", format='jpg')
      plt.show()
    
    
    
    def plot_confusion_matrix(self,save=False):
        
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
            
        cm = confusion_matrix(self.y, self.ypred)
        confusion = sns.heatmap(cm, annot = True,fmt='.4g')
        if save is True:
            # Save the pairplot as jpg
            confusion.savefig("confusion_matrix.jpg", format='jpg')

        # Show the pairplot
        plt.show()

         
    '''
   
        
      y_pred = gs_svc.predict(X).astype(int)
      y_score = gs_svc.predict_proba(X)
    
    plot_roc_and_precision_recall(Y, y_score[:,1])  # provide the column for the scores belonging only to the positive class

    '''

    def plot_precision_recall_cv(self,save=False):
        #some chekc  for models predict proba AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
        
        
        y_scores = cross_val_predict(self.best_estimator , self.X, self.y, cv=5, method='predict_proba')
        precision, recall, _ = precision_recall_curve(self.y, y_scores[:, 1])
        average_precision = average_precision_score(self.y, y_scores[:, 1])

        plt.figure()
        plt.step(recall, precision, where='post', color='b', alpha=0.5)
        plt.fill_between(recall, precision, step='post', alpha=0.5, color='b')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("precision_recall_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()


    def plot_roc_curve_cv(self,save=False):
        #also some chekc for some models AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
            
        y_scores = cross_val_predict(self.best_estimator, self.X, self.y, cv=5, method='predict_proba')
        fpr, tpr, _ = roc_curve(self.y, y_scores[:, 1])
        roc_auc = auc(fpr, tpr)

        plt.figure()
        lw = 2
        plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend(loc="lower right")
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("roc_curve_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()



    def intercept(self):
        if self.best_estimator.intercept_ is None or self.best_estimator.coef_ is None :
            raise AttributeError('Your model can not provide coefficients and intercept, Make sure you use appropriate model!')

        print("intercept_:")
        print(self.best_estimator.intercept_ )
        print()
        print("coef_:")
        print(self.best_estimator.coef_)
        
        
        
    def plot_feature_importance(self,save=False):
        
        if self.best_estimator.feature_importances_ is None:
            raise AttributeError('This function is used only for Random Forest classifier!')
            
            
        dataset=load_breast_cancer()
        n_feature = 30
        plt.barh(range(n_feature), self.best_estimator.feature_importances_, align='center')
        plt.yticks(np.arange(n_feature), dataset.feature_names)
        plt.xlabel('Feature Importance')
        plt.ylabel('Feature')
        if save is True:
            plt.savefig("feature_importance.jpg", format='jpg')
        plt.show()
        
        
    
    def plot_coeff_weitgh(self,save=False):
        if self.best_estimator.coefs_ is None:
            raise AttributeError('This function is used only for MLP classifier!')
            
        dataset=load_breast_cancer()
        
        plt.figure(figsize=(20,5))
        plt.imshow(self.best_estimator.coefs_[0],interpolation='None',cmap='GnBu')
        plt.yticks(range(30), dataset.feature_names)
        plt.xlabel('Colums in weight matrix')
        plt.ylabel('Input feature')
        plt.colorbar()
        if save is True:
            plt.savefig("coeff_weight.jpg", format='jpg')
        plt.show()
        


    def plot_true_vs_predict(self,name_model,save=False):
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
        # Perform PCA for dimensionality reduction
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(self.X)

        # Perform cross-validation prediction
        y_pred = cross_val_predict(self.best_estimator, self.X, self.y, cv=5)
        
        # Visualize the true labels versus predicted labels for training set
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=np.array(self.y), cmap='viridis', label='True Labels')
        plt.title(f'True Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=self.ypred, cmap='viridis', label='Predicted Labels')
        plt.title(f'Predicted Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.tight_layout()
        if save is True:
            plt.savefig("true_vs_predict.jpg", format='jpg')
        plt.show()
        


        
    def plot_permutation_importance(self,save=False):
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
        
        # Calculate Permutation Feature Importance
        result = permutation_importance(self.best_estimator, self.X, self.y, n_repeats=10, random_state=42, n_jobs=2)
        perm_sorted_idx = result.importances_mean.argsort()
        
        fig, ax = plt.subplots(figsize=(8, 6))  # Create a new figure and axes
        ax.boxplot(
            result.importances[perm_sorted_idx].T,
            vert=False,
            labels=pd.DataFrame(self.X).columns[perm_sorted_idx],
            )
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Permutation Importance")
        ax.set_title("Permutation Feature Importance")
        if save is True:
            plt.savefig("permutation_importance.jpg", format='jpg')
        plt.show()











class classification_PIPE:

    def __init__(self, X, y):
        
        if X is None or y is None:
            raise ValueError("Insert Both x and y as input")
            
        self.X=X
        self.y= np.array(y).ravel()
        
        
        self.model=None
        self.ypred=None
        self.best_estimator=None
        self.gs=None
        self.myparams=None
        self.pipe=None
        
        
        #self.best_estimator.named_steps['model'].intercept_ = None  
        #self.best_estimator.named_steps['model'].coef_ =None  
        #self.best_estimator.named_steps['model'].feature_importances_ = None
        
        
        
    def data_clear_check(self):
       self.df=pd.concat([pd.DataFrame(self.X),self.y],axis=1)
       
       #df=self.df
       print('-------------------')
       self.df.isnull().sum()
       print('-------------------')
       self.df.isna().sum()
       print('-------------------')
       self.df.describe()
       
       
    def plot_correlation(self,save=False):
        yX_corr = self.df.corr(method='pearson')

          # Convert to abolute values
        yX_abs_corr = np.abs(yX_corr) 
          

        plt.figure(figsize=(10, 10))
        plt.imshow(yX_abs_corr, cmap='RdYlGn', interpolation='none', aspect='auto')
        plt.colorbar()
        plt.xticks(range(len(yX_abs_corr)), yX_abs_corr.columns, rotation='vertical')
        plt.yticks(range(len(yX_abs_corr)), yX_abs_corr.columns);
        plt.suptitle('Pearson Correlation Heat Map (absolute values)', fontsize=15, fontweight='bold')
        if save is True:
            plt.savefig("correlation.jpg", format='jpg')
        plt.show()
        
            
    def plot_pair_plot(self,save=False):
        pairplot = sns.pairplot(self.df, hue='diognosis', palette='coolwarm')
        if save is True:
            # Save the pairplot as jpg
            pairplot.savefig("pairplot.jpg", format='jpg')

        # Show the pairplot
        plt.show()

    def make_pipe(self,model,poly=None,specific_preprocessor=None):
        if poly!=None and specific_preprocessor==None:
            self.model=model
            self.scaler=MinMaxScaler()
            self.poly=PolynomialFeatures()
            self.pipe=Pipeline([("scaler", self.scaler), ('poly',self.poly) ,("model", self.model)])
            
            print('Pipeline is constructed')

        elif specific_preprocessor!=None and poly==None:
            self.model=model
            self.scaler=specific_preprocessor

            self.pipe=Pipeline([("scaler", self.scaler) ,("model", self.model)])
            print('Pipeline is constructed')
            

         
        elif specific_preprocessor!=None and poly!=None:
            self.model=model
            self.scaler=specific_preprocessor
            self.poly=PolynomialFeatures()

            self.pipe=Pipeline([("scaler", self.scaler), ('poly',self.poly) ,("model", self.model)])
            print('Pipeline is constructed')
            
            
        
        else:
            self.model=model
            self.scaler=MinMaxScaler()
            self.pipe=Pipeline([("scaler", self.scaler), ("model", self.model)])
            
            print('Pipeline is constructed')
          
         
         
    def fitting(self,myparams):
        if self.pipe is None:
            raise ValueError('please first use make_pipe fucntion!')
            

        self.myparams=myparams
        
        self.gs= GridSearchCV(self.pipe, self.myparams , return_train_score=True)
        self.gs.fit(self.X,self.y)
        
        self.best_estimator=self.gs.best_estimator_
        
        self.ypred=cross_val_predict(self.best_estimator,self.X,self.y)
        
        print('The model is trained and is ready now!')
        
        
    def best(self):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
            
        self.best_estimator=self.gs.best_estimator_
        self.best_score=self.gs.best_score_
        return self.best_estimator,self.best_score
    
    
    def predicting(self,new_x):
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
                
        prediction=self.best_estimator.predict(new_x)
        return prediction

    def plot_grid_search(self, save=False):
      '''
      Plot 2d grid search heatmap
      
      Parameters:
      grid_search: instance of sklearn's GridSearchCV
      grid_params: dictionary of grid search parameters  
      x_param: name of x-axis parameter in grid_params
      y_param: name of y-axis parameter in grid_params
      is_verbose (optional): print results
      
      Return:
      grid_search.best_score_: best score found
      grid_search.best_estimator_: best estimator found
      '''
      
      #check because we have scaler here

      key_params=list(self.myparams.keys())
      
      grid_params_x = self.myparams[key_params[0]]
      grid_params_y = self.myparams[key_params[1]]
      
      if self.gs is None:
          raise ValueError("Model not trained yet. Call fitting() first.")
          
      if self.myparams is None:
          raise ValueError("Model not trained yet. Call fitting() first.")
      
        
      df_results = pd.DataFrame(self.gs.cv_results_)
      ar_scores = np.array(df_results.mean_test_score).reshape(len(grid_params_y), len(grid_params_x))
      sns.heatmap(ar_scores, annot=True, fmt='.3f', xticklabels=grid_params_x, yticklabels=grid_params_y)
      print()
      plt.suptitle('Grid search heatmap')
      plt.xlabel(key_params[0])
      plt.ylabel(key_params[1])
      if save is True:
          # Save the pairplot as jpg
          plt.savefig("grid_search_heat_map.jpg", format='jpg')
      plt.show()
    
    
    
    def plot_confusion_matrix(self,save=False):
        
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
            
        cm = confusion_matrix(self.y, self.ypred)
        confusion = sns.heatmap(cm, annot = True,fmt='.4g')
        if save is True:
            # Save the pairplot as jpg
            confusion.savefig("confusion_matrix.jpg", format='jpg')

        # Show the pairplot
        plt.show()

         
    '''
   
        
      y_pred = gs_svc.predict(X).astype(int)
      y_score = gs_svc.predict_proba(X)
    
    plot_roc_and_precision_recall(Y, y_score[:,1])  # provide the column for the scores belonging only to the positive class

    '''

    def plot_precision_recall_cv(self,save=False):
        #some chekc  for models predict proba AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
        
        
        y_scores = cross_val_predict(self.best_estimator , self.X, self.y, cv=5, method='predict_proba')
        precision, recall, _ = precision_recall_curve(self.y, y_scores[:, 1])
        average_precision = average_precision_score(self.y, y_scores[:, 1])

        plt.figure()
        plt.step(recall, precision, where='post', color='b', alpha=0.5)
        plt.fill_between(recall, precision, step='post', alpha=0.5, color='b')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("precision_recall_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()


    def plot_roc_curve_cv(self,save=False):
        #also some chekc for some models AttributeError:
        if self.model is None:
            raise ValueError("Model not trained yet. Call fitting() first.")
            
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
            
        y_scores = cross_val_predict(self.best_estimator, self.X, self.y, cv=5, method='predict_proba')
        fpr, tpr, _ = roc_curve(self.y, y_scores[:, 1])
        roc_auc = auc(fpr, tpr)

        plt.figure()
        lw = 2
        plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend(loc="lower right")
        if save is True:
            # Save the pairplot as jpg
            plt.savefig("roc_curve_cv.jpg", format='jpg')

        # Show the pairplot
        plt.show()



    def intercept(self):
        if self.best_estimator.named_steps['model'].intercept_ is None or self.best_estimator.named_steps['model'].coef_ is None :
            raise AttributeError('Your model can not provide coefficients and intercept, Make sure you use appropriate model!')

        print("intercept_:")
        print(self.best_estimator.named_steps['model'].intercept_ )
        print()
        print("coef_:")
        print(self.best_estimator.named_steps['model'].coef_)
        
        
        
    def plot_feature_importance(self,save=False):
        
        if self.best_estimator.named_steps['model'].feature_importances_ is None:
            raise AttributeError('This function is used only for Random Forest classifier!')
            
            
        dataset=load_breast_cancer()
        n_feature = 30
        plt.barh(range(n_feature), self.best_estimator.named_steps['model'].feature_importances_, align='center')
        plt.yticks(np.arange(n_feature), dataset.feature_names)
        plt.xlabel('Feature Importance')
        plt.ylabel('Feature')
        if save is True:
            plt.savefig("feature_importance.jpg", format='jpg')
        plt.show()
        
        
    
    def plot_coeff_weitgh(self,save=False):
        if self.best_estimator.named_steps['model'].coefs_ is None:
            raise AttributeError('This function is used only for MLP classifier!')
            
        dataset=load_breast_cancer()
        
        plt.figure(figsize=(20,5))
        plt.imshow(self.best_estimator.named_steps['model'].coefs_[0],interpolation='None',cmap='GnBu')
        plt.yticks(range(30), dataset.feature_names)
        plt.xlabel('Colums in weight matrix')
        plt.ylabel('Input feature')
        plt.colorbar()
        if save is True:
            plt.savefig("coeff_weight.jpg", format='jpg')
        plt.show()
        


    def plot_true_vs_predict(self,name_model,save=False):
        if self.ypred is None:
            raise ValueError('Predicted value is not fined, Call Predicting first.')
            
        # Perform PCA for dimensionality reduction
        #new pipe
        
        #self.pca = PCA(n_components=2)
        
        #self.pca_pipe=Pipeline([("scaler", self.scaler), ('pca',self.pca),("model", self.model)])
        

        # Perform cross-validation prediction
        #self.y_pca_pred = cross_val_predict(self.pca_pipe, self.X, self.y, cv=5)
        
        
        pca = PCA(n_components=2)
        scaler=MinMaxScaler()

        xstandard=scaler.fit.transform(self.X)
        
        X_pca = pca.fit_transform(xstandard)
        
        
        
        # Visualize the true labels versus predicted labels for training set
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=np.array(self.y), cmap='viridis', label='True Labels')
        plt.title(f'True Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=self.y_pred, cmap='viridis', label='Predicted Labels')
        plt.title(f'Predicted Labels for model: {name_model}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()

        plt.tight_layout()
        if save is True:
            plt.savefig("true_vs_predict.jpg", format='jpg')
        plt.show()
        


        
    def plot_permutation_importance(self,save=False):
        if self.best_estimator is None:
            raise ValueError("Best is not fined. Call best() first.")
        
        # Calculate Permutation Feature Importance
        result = permutation_importance(self.best_estimator, self.X, self.y, n_repeats=10, random_state=42, n_jobs=2)
        perm_sorted_idx = result.importances_mean.argsort()
        
        fig, ax = plt.subplots(figsize=(8, 6))  # Create a new figure and axes
        ax.boxplot(
            result.importances[perm_sorted_idx].T,
            vert=False,
            labels=pd.DataFrame(self.X).columns[perm_sorted_idx],
            )
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Permutation Importance")
        ax.set_title("Permutation Feature Importance")
        if save is True:
            plt.savefig("permutation_importance.jpg", format='jpg')
        plt.show()





