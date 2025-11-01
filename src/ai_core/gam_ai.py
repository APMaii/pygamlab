import json, base64, joblib, io, os

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    confusion_matrix, ConfusionMatrixDisplay,
    classification_report, r2_score, mean_absolute_error, mean_squared_error
)


  

class Gam_Ai_Workflow:
    """
    Intelligent, type-aware workflow manager for `.gam_ai` model packages.

    This class provides a unified interface for loading, evaluating, refitting,
    and visualizing machine learning models saved in the `.gam_ai` format.
    It automatically detects the model type (classifier, regressor, or
    unsupervised) and routes the evaluation pipeline accordingly.

    Parameters
    ----------
    model_name : str
        The name (without extension) of the `.gam_ai` file to be loaded.
    base_dir : str, optional, default="gam_models"
        Directory containing saved `.gam_ai` model files.

    Raises
    ------
    FileNotFoundError
        If the specified `.gam_ai` file cannot be located within `base_dir`.

    Notes
    -----
    This class depends on the `GamAI_io` handler for deserializing `.gam_ai` files.
    Each `.gam_ai` file contains both model metadata and a serialized model object.
    Once loaded, `Gam_Ai_Workflow` provides:
        - Smart evaluation (`evaluate()`)
        - Visualization (`visualize_unsupervised()`)
        - Refit capabilities (`refit()`)
        - Summaries (`summary()`)

    Examples
    --------
    >>> workflow = Gam_Ai_Workflow("rf_classifier_v1", base_dir="models")
    ✅ Loaded model 'rf_classifier_v1' (classifier) successfully.
    >>>
    >>> workflow.summary()
    📘 MODEL SUMMARY
    model_name: rf_classifier_v1
    model_type: classifier
    author_name: John Doe
    best_accuracy: 0.94
    ...
    >>>
    >>> workflow.evaluate(X_test, y_test)
    🎯 Accuracy: 0.9470
    📊 Classification Report:
    ...
    """
    def __init__(self, model_name, base_dir=None):
        self.model_name = model_name


        self.gam=GAM_AI_MODEL(model_name,base_dir=base_dir)
        self.ml_model = self.gam.ml_model
        self.model_type = self.gam.model_type.lower()

        print(f"✅ Loaded model '{self.model_name}' ({self.model_type}) successfully.")

    # ---------- General ----------
    def summary(self):
        """
        Display a concise summary of the loaded model and its metadata.

        Prints all metadata fields stored in the `.gam_ai` file,
        including model type, author information, and hyperparameters.

        Returns
        -------
        None
        """
        print("\n📘 MODEL SUMMARY")
        for k, v in self.gam.__dict__.items():
            if k != "model":
                print(f"{k}: {v}")

    def predict(self, X):
        """
        Run model inference on input data.

        Parameters
        ----------
        X : array-like
            Input features compatible with the trained model.

        Returns
        -------
        np.ndarray
            Model predictions corresponding to `X`.

        Raises
        ------
        ValueError
            If no model is loaded.
        """
        if self.ml_model is None:
            raise ValueError("Model not loaded.")
        return self.ml_model.predict(X)

    def refit(self, X, y):
        """
        Retrain (refit) the loaded model on new data.

        Parameters
        ----------
        X : array-like
            Training features.
        y : array-like
            Corresponding training labels or targets.

        Raises
        ------
        NotImplementedError
            If the model does not support the `.fit()` method.

        Notes
        -----
        This method modifies the current model in place and does not automatically
        update the `.gam_ai` file on disk. To persist changes, re-save the model
        using `GamAI_io.save()` after refitting.
        """
        if hasattr(self.ml_model, "fit"):
            self.ml_model.fit(X, y)
            print("🔁 Model refitted successfully.")
        else:
            raise NotImplementedError("This model type cannot be refitted.")

    # ---------- Classifier ----------
    def evaluate_classifier(self, X, y_true):
        """
        Evaluate a classification model and display key performance metrics.

        Parameters
        ----------
        X : array-like
            Input test features.
        y_true : array-like
            Ground-truth class labels.

        Returns
        -------
        None

        Notes
        -----
        This method computes and displays:
            - Accuracy score
            - Classification report
            - Confusion matrix plot
        """
        
        y_pred = self.ml_model.predict(X)
        acc = np.mean(y_pred == y_true)
        print(f"🎯 Accuracy: {acc:.4f}\n")
        print("📊 Classification Report:\n")
        print(classification_report(y_true, y_pred))

        cm = confusion_matrix(y_true, y_pred)
        disp = ConfusionMatrixDisplay(cm)
        disp.plot(cmap="Blues")
        plt.title(f"Confusion Matrix: {self.gam.model_name}")
        plt.show()

    # ---------- Regressor ----------
    def evaluate_regressor(self, X, y_true):
        """
        Evaluate a regression model and visualize performance trends.

        Parameters
        ----------
        X : array-like
            Input test features.
        y_true : array-like
            True continuous target values.

        Returns
        -------
        None

        Notes
        -----
        This method prints and plots:
            - Coefficient of determination (R²)
            - Mean Absolute Error (MAE)
            - Mean Squared Error (MSE)
            - Predicted vs. True value scatter plot
        """
        y_pred = self.ml_model.predict(X)
        r2 = r2_score(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        mse = mean_squared_error(y_true, y_pred)

        print(f"📈 R²: {r2:.4f}")
        print(f"📉 MAE: {mae:.4f}")
        print(f"📉 MSE: {mse:.4f}\n")

        plt.figure(figsize=(6, 4))
        plt.scatter(y_true, y_pred, alpha=0.6)
        plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'r--')
        plt.xlabel("True Values")
        plt.ylabel("Predicted Values")
        plt.title(f"Prediction Trend — {self.gam.model_name}")
        plt.grid(True)
        plt.show()

    # ---------- Unsupervised ----------
    def visualize_unsupervised(self, X):
        """
        Visualize cluster assignments or feature transformations for unsupervised models.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data to visualize.

        Returns
        -------
        None

        Raises
        ------
        NotImplementedError
            If the model lacks both `predict()` and `transform()` methods.

        Notes
        -----
        - If the model has a `.predict()` method, cluster assignments are plotted.
        - If the model has a `.transform()` method, the transformed feature space is shown.
        - This visualization assumes that the first two components or features
          are suitable for 2D projection.
        """
        if hasattr(self.ml_model, "predict"):
            y_pred = self.ml_model.predict(X)
            plt.scatter(X[:, 0], X[:, 1], c=y_pred, cmap="viridis", s=30)
            plt.title(f"Cluster Visualization — {self.gam.model_name}")
            plt.show()
        elif hasattr(self.ml_model, "transform"):
            X_trans = self.ml_model.transform(X)
            plt.scatter(X_trans[:, 0], X_trans[:, 1], s=30)
            plt.title(f"Feature Space — {self.gam.model_name}")
            plt.show()
        else:
            raise NotImplementedError("This unsupervised model has no visualization method.")

    # ---------- Smart Evaluation ----------
    def evaluate(self, X, y_true=None):
        """
        Automatically dispatch model evaluation based on its declared type.

        Parameters
        ----------
        X : array-like
            Input data.
        y_true : array-like, optional
            Ground-truth labels or target values (required for supervised models).

        Raises
        ------
        ValueError
            If the model type is not recognized.

        Notes
        -----
        This method intelligently determines which evaluation routine to run:
            - `evaluate_classifier()` for classification models
            - `evaluate_regressor()` for regression models
            - `visualize_unsupervised()` for unsupervised models

        Examples
        --------
        >>> workflow.evaluate(X_test, y_test)
        🎯 Accuracy: 0.9470
        📊 Classification Report:
        ...
        """
        if self.model_type == "classifier":
            self.evaluate_classifier(X, y_true)
        elif self.model_type == "regressor":
            self.evaluate_regressor(X, y_true)
        elif self.model_type == "unsupervised":
            self.visualize_unsupervised(X)
        else:
            raise ValueError(f"Unknown model_type: {self.model_type}")



    def get_GAM_AI_MODEL(self):
        """
        Retrieve the fully loaded `.gam_ai` model instance associated with this workflow.

        Returns
        -------
        GAM_AI_MODEL
            The `GAM_AI_MODEL` object currently managed by this workflow instance.  
            This object encapsulates both the model’s metadata (e.g., author info, 
            training details, performance metrics) and the deserialized scikit-learn model 
            accessible via the attribute `ml_model`.

        Examples
        --------
        >>> workflow = Gam_Ai_Workflow("cu-nanocomposites-poisson-ratio-lr")
        ✅ Loaded model 'cu-nanocomposites-poisson-ratio-lr' (train/test) successfully.

        >>> gam_model = workflow.get_GAM_AI_MODEL()
        >>> type(gam_model)
        <class 'PyGamLab.ai_core.gam_ai.GAM_AI_MODEL'>

        >>> gam_model.summary()
        📘 MODEL METADATA SUMMARY
        model_name: cu-nanocomposites-poisson-ratio-lr
        author_name: Shaoyu Zhao, Yingyan Zhang, Yihe Zhang et al.
        best_accuracy: {'MAE': 0.0541, 'MSE': 0.0042, 'R2': 0.39}
        ⚙️ ML Model: <class 'sklearn.linear_model._base.LinearRegression'>

        Notes
        -----
        This method serves as a safe accessor for the underlying `GAM_AI_MODEL` instance 
        (`self.gam`) loaded during `Gam_Ai_Workflow` initialization.  
        It can be used to directly inspect model metadata, retrieve the raw ML model, 
        or perform low-level analysis without invoking higher-level workflow methods.
        """
        return self.gam


            
            







import os, json, io, base64, joblib

class GAM_AI_MODEL:
    """
    A unified data structure for loading and representing `.gam_ai` model packages.

    The `.gam_ai` format encapsulates both the machine learning model (serialized
    via `joblib` and encoded in Base64) and its accompanying metadata. 
    This class provides a standardized interface for deserializing, inspecting,
    and utilizing such packaged models within the PyGamLab AI ecosystem.

    The class dynamically attaches metadata fields (e.g., `author_name`, 
    `model_type`, `best_accuracy`, etc.) as instance attributes and exposes 
    the trained ML model under the attribute `ml_model`.

    Parameters
    ----------
    model_name : str
        The name of the model file (without the `.gam_ai` extension) to be loaded.

    base_dir : str, optional
        The base directory containing `.gam_ai` model files. If not provided,
        the class automatically defaults to the `gam_models` directory located
        alongside this module file.

    Attributes
    ----------
    model_name : str
        The identifier of the loaded `.gam_ai` model.
    
    file_path : str
        Full path to the `.gam_ai` file on disk.
    
    ml_model : object
        The deserialized scikit-learn model instance (e.g., LinearRegression, RandomForestRegressor, etc.).
    
    <dynamic_metadata_fields> : Any
        All key-value pairs from the `"metadata"` section of the `.gam_ai` file are 
        dynamically added as attributes (e.g., `author_name`, `description`, `best_accuracy`, etc.).

    Raises
    ------
    FileNotFoundError
        If the specified `.gam_ai` file cannot be found at the resolved `file_path`.

    Examples
    --------
    >>> from PyGamLab.ai_core.gam_ai import GAM_AI_MODEL
    >>> model = GAM_AI_MODEL("cu-nanocomposites-poisson-ratio-lr")
    ✅ Loaded GAM_AI_MODEL: 'cu-nanocomposites-poisson-ratio-lr'

    >>> model.summary()
    📘 MODEL METADATA SUMMARY
    model_name: cu-nanocomposites-poisson-ratio-lr
    model_type: train/test
    author_name: Shaoyu Zhao, Yingyan Zhang, Yihe Zhang et al.
    best_accuracy: {'MAE': 0.0541, 'MSE': 0.0042, ...}
    ⚙️ ML Model: <class 'sklearn.linear_model._base.LinearRegression'>

    Notes
    -----
    - The `.gam_ai` file format is designed to preserve reproducibility of ML experiments.
    - Metadata is stored in JSON format, while the ML model itself is serialized using `joblib`
      and encoded with Base64 for portability.
    - This class does not train or evaluate models; it only loads and interprets
      pre-trained model artifacts.

    See Also
    --------
    Gam_Ai_Workflow : High-level workflow manager that builds upon this class to 
                      provide evaluation, visualization, and retraining utilities.
    """

    def __init__(self, model_name, base_dir=None):
        self.model_name = model_name

        # Default directory if none provided
        if base_dir is None:
            base_dir = os.path.join(os.path.dirname(__file__), "gam_models")

        self.file_path = os.path.join(base_dir, f"{model_name}.gam_ai")

        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"❌ File not found: {self.file_path}")

        # --- Load JSON file ---
        with open(self.file_path, "r") as f:
            data = json.load(f)

        metadata = data.get("metadata", {})
        model_data = data.get("model_data", None)

        # --- Assign metadata attributes dynamically ---
        for key, value in metadata.items():
            setattr(self, key, value)

        # --- Decode Base64 model ---
        if model_data:
            model_bytes = base64.b64decode(model_data)
            buffer = io.BytesIO(model_bytes)
            self.ml_model = joblib.load(buffer)
        else:
            self.ml_model = None

        print(f"✅ Loaded GAM_AI_MODEL: '{self.model_name}'")

    def summary(self):
        """
        Print a structured overview of all model metadata.

        This method displays the dynamically loaded metadata attributes 
        (e.g., author information, accuracy metrics, DOI, etc.) and the 
        associated machine learning model type.

        Examples
        --------
        >>> model.summary()
        📘 MODEL METADATA SUMMARY
        model_name: cu-nanocomposites-poisson-ratio-lr
        model_type: train/test
        author_name: Shaoyu Zhao
        best_accuracy: {'MAE': 0.054, 'MSE': 0.0042}
        ⚙️ ML Model: <class 'sklearn.linear_model._base.LinearRegression'>
        """

        print("\n📘 MODEL METADATA SUMMARY\n")
        for k, v in self.__dict__.items():
            if k not in ("ml_model", "file_path"):
                print(f"{k}: {v}")
        print("\n⚙️ ML Model:", type(self.ml_model))






class GamAI_io:
    """
    A unified input/output handler for saving and loading `.gam_ai` model packages.

    This class encapsulates both machine learning model serialization and
    relevant metadata (e.g., author info, model parameters, training details)
    into a single portable `.gam_ai` file. It enables seamless model deployment,
    archival, and reproducibility by combining the binary model object and
    human-readable metadata in one structured JSON container.

    Parameters
    ----------
    model_name : str, optional
        A short, descriptive name for the model. Used as the filename during saving.
    model_type : str, optional
        The algorithmic or architectural family of the model
        (e.g., "RandomForest", "XGBoost", "CNN", "Transformer").
    description : str, optional, default=""
        A brief summary describing the model’s purpose, training dataset,
        or key features.
    author_name : str, optional, default=""
        Full name of the model creator.
    author_email : str, optional, default=""
        Contact email for correspondence or citation.
    trainer_name : str, optional, default=""
        The individual or system responsible for model training.
    best_accuracy : float, optional
        The highest validation or test accuracy achieved during training.
    doi : str, optional
        Digital Object Identifier (DOI) associated with the published model or dataset.
    hyperparam_range : dict, optional, default={}
        Dictionary defining hyperparameter search ranges used during optimization.
    best_params : dict, optional, default={}
        Dictionary of the final optimized hyperparameter values.
    ml_model : object, optional
        The trained machine learning model instance (e.g., scikit-learn estimator).

    Notes
    -----
    The `.gam_ai` file format consists of a JSON object with two main sections:
    
    - **metadata**: Contains descriptive fields such as author, model type,
      and hyperparameters.
    - **model_data**: Contains the Base64-encoded binary serialization of the
      trained model, produced via `joblib`.

    This approach ensures full portability and JSON readability, enabling
    both programmatic and manual inspection of model metadata.

    Examples
    --------
    >>> from gamai_io import GamAI_io
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> model = RandomForestClassifier(n_estimators=100, random_state=42)
    >>> model.fit(X_train, y_train)
    >>> 
    >>> package = GamAI_io(
    ...     model_name="rf_classifier_v1",
    ...     model_type="RandomForest",
    ...     description="Predicts material phases using compositional data",
    ...     author_name="John Doe",
    ...     author_email="john.doe@example.com",
    ...     best_accuracy=0.94,
    ...     hyperparam_range={"n_estimators": [50, 100, 200]},
    ...     best_params={"n_estimators": 100, "max_depth": 10},
    ...     ml_model=model
    ... )
    >>> 
    >>> # Save to file
    >>> package.save(save_dir="models")
    💾 Saved ml_model package: models/rf_classifier_v1.gam_ai
    >>>
    >>> # Load from file
    >>> loaded_package = GamAI_io.load("models/rf_classifier_v1.gam_ai")
    >>> restored_model = loaded_package.ml_model
    >>> restored_model.predict(X_test[:5])

    See Also
    --------
    joblib.dump : Efficient serialization of Python objects.
    json : Standard JSON encoder/decoder.
    base64 : Encoding binary model data for safe JSON storage.

    """

    def __init__(self, **kwargs):
        self.model_name = kwargs.get("model_name")
        self.model_type = kwargs.get("model_type")
        self.description = kwargs.get("description", "")
        self.author_name = kwargs.get("author_name", "")
        self.author_email = kwargs.get("author_email", "")
        self.trainer_name = kwargs.get("trainer_name", "")
        self.best_accuracy = kwargs.get("best_accuracy", None)
        self.doi=kwargs.get("doi", None)
        self.hyperparam_range = kwargs.get("hyperparam_range", {})
        self.best_params = kwargs.get("best_params", {})
        self.ml_model = kwargs.get("ml_model", None)

    # ----------- Save ----------
    def save(self, save_dir="models"):
        """
        Serialize and save the current model and metadata as a `.gam_ai` package.

        Parameters
        ----------
        save_dir : str, optional, default="models"
            The target directory to save the `.gam_ai` file.
            The directory will be created if it does not exist.

        Raises
        ------
        ValueError
            If no model (`ml_model`) is attached to the current instance.

        Notes
        -----
        The model is serialized via `joblib` and encoded with Base64 to
        ensure JSON compatibility. The resulting file can be safely shared
        or uploaded to repositories without binary corruption.

        """
        os.makedirs(save_dir, exist_ok=True)
        filename = f"{self.model_name}.gam_ai"
        filepath = os.path.join(save_dir, filename)

        if self.ml_model is None:
            raise ValueError("No ml_model attached to save inside .gam_ai.")

        buffer = io.BytesIO()
        joblib.dump(self.ml_model, buffer)
        buffer.seek(0)
        encoded_model = base64.b64encode(buffer.read()).decode("utf-8")

        data = {
            "metadata": {k: v for k, v in self.__dict__.items() if k != "ml_model"},
            "model_data": encoded_model
        }

        with open(filepath, "w") as f:
            json.dump(data, f, indent=4)
        print(f"💾 Saved ml_model package: {filepath}")

    # ----------- Load ----------
    @staticmethod
    def load(filepath):
        """
        Load a `.gam_ai` model package from disk.

        Parameters
        ----------
        filepath : str
            Full path to the `.gam_ai` file to be loaded.

        Returns
        -------
        GamAI_io
            A `GamAI_io` instance containing both metadata and the
            deserialized machine learning model (`ml_model`).

        Notes
        -----
        The loading process reverses the Base64 encoding and `joblib`
        serialization to reconstruct the original model object.

        """
        with open(filepath, "r") as f:
            data = json.load(f)
        metadata = data["metadata"]
        model_data = base64.b64decode(data["model_data"])
        buffer = io.BytesIO(model_data)
        model = joblib.load(buffer)
        return GamAI_io(**metadata, model=model)
    
    








