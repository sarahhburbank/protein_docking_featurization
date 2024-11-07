import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
import pandas as pd
import ast
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

def preprocess(results_path):
    #header = ['name', 'crmsd', 'Predicted Hinge Count',  'kingdom', 'Synthetic', 'Number of Ligand Atoms','Number of Receptor Atoms','   Ratio of(# Receptor Atoms):(# of Ligand Atoms)', 'Number of Ligand Chains', 'Number of Receptor Chains', 'Ratio of(# Receptor Chains):(# of Ligand Chains)', 'Std Dev Ligand Atoms/Chain', 'the Std Dev Receptor Atoms/Chain', 'Ratio of Disulfide Bonds to Receptor Atoms']
    header = ['name', 'crmsd', 'Hinge Count',  'kingdom', 'Synthetic', '# Lg. Atoms','# Rc. Atoms','Ratio Atoms', ' # Lg.  Chains', ' # Rc. Chains', 'Chains Ratio', 'Std Dev L.', 'Std Dev R.', 'DiSulf. Ratio']

    # Load the CSV file
    data = pd.read_csv(results_path)
    # Print the number of columns

    data.columns = header
    for index, __ in data.iterrows():
        if data.at[index, "crmsd"] < 9.41:
            
            data.at[index, "crmsd"] = 1
        elif data.at[index, "crmsd"] >= 9.41 and data.at[index, "crmsd"] < 21.03:
            
            data.at[index, "crmsd"] = 2
        elif data.at[index, "crmsd"] >= 21: #and data.at[index, "crmsd"] < 25:
            data.at[index, "crmsd"] = 3
        ''' elif data.at[index, "crmsd"] >= 25:
            data.at[index, "crmsd"] = 4'''
       
            
    print(type(data['kingdom']))
    data['kingdom'] = data['kingdom'].apply(lambda x: list(map(int, x.strip('()').split())))
    data_2d = data['kingdom'].apply(pd.Series) #shape is ask expceted 344,4
    new_columns = ['Eukaryotic', 'Prokaryotic', 'Viral', 'Other']

# Change the header
    data_2d.columns = new_columns
    data_encoded = pd.concat([data.drop('kingdom', axis=1), data_2d], axis=1)

    exclude_columns = ['name', 'crmsd']

    # Get columns to normalize
    columns_to_normalize = [col for col in data_encoded.columns if col not in exclude_columns]

    # Normalize columns except for 'A' and 'B'
    scaler = MinMaxScaler(feature_range=(-1, 1))

    data_encoded.columns = data_encoded.columns.astype(str)
    print(data_encoded.columns)
    print("columns to normalize", columns_to_normalize)
    data_encoded[columns_to_normalize] = scaler.fit_transform(data_encoded[columns_to_normalize])
    print(data_encoded)
    return data_encoded
import wandb
def rfr_classifier(data):
#['crmsd', 'name',  'Hinge Count', 'DiSulf. Ratio', 'Ratio Atoms', 'Chains Ratio'
    X = data.drop(columns=['crmsd', 'name',  'Eukaryotic', 'Prokarytoic', 'Viral', 'Other'])
    y = data['crmsd']
    feature_labels = X.columns  # Get the feature labels
    print(feature_labels)

    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    train_accuracy_history = []
    test_accuracy_history = []
   # Initialize the Random Forest Classifier
    for i in range(1, 100):  # Train for 100 iterations
        rfc = RandomForestClassifier(n_estimators=i, max_depth=13, min_samples_split=5, min_samples_leaf=5,)
        rfc.fit(X_train, y_train)

        # Calculate train and test accuracy
        train_pred = rfc.predict(X_train)
        test_pred = rfc.predict(X_test)
        train_accuracy = accuracy_score(y_train, train_pred)
        test_accuracy = accuracy_score(y_test, test_pred)

        # Append accuracy scores to history lists
        train_accuracy_history.append(train_accuracy)
        test_accuracy_history.append(test_accuracy)

# Plot the change in test and train accuracy as the model runs
    
    print("train", train_accuracy)
    print("test", test_accuracy)
    from scipy.interpolate import make_interp_spline

# Plot the change in test and train accuracy as the model runs
    # Create smoother lines by interpolating the data
    #x_smooth = np.linspace(1, 200, 400)  # Increase the number of points for smoother curves
    #train_smooth = make_interp_spline(range(1, 100), train_accuracy_history)(x_smooth)
    #test_smooth = make_interp_spline(range(1, 100), test_accuracy_history)(x_smooth)

    # Plot the smoothed lines
    plt.figure(figsize=(10, 6))

    # Plot smoothed train accuracy
    #plt.plot(x_smooth, train_smooth, label='Train Accuracy', color='#1f77b4')  # Blue color

    # Plot smoothed test accuracy
   # plt.plot(x_smooth, test_smooth, label='Test Accuracy', color='#ff7f0e')  # Orange color

    # Add labels and title
   # plt.xlabel('Number of Trees')
    #plt.ylabel('Accuracy')
    #plt.title('Num Trees vs Test and Train Accuracy ')

    # Add legend
   #plt.legend()

    # Show the plot
   # plt.show()

    # Display feature importance
    feature_importance = pd.Series(rfc.feature_importances_, index=feature_labels)  # Use feature labels as index
    feature_importance_sorted = feature_importance.sort_values(ascending=False) 
    feature_importance_sorted.plot(kind='bar')

    plt.xlabel('Feature',fontsize=6)
    plt.ylabel('Importance')
    plt.title('Feature Importance')
    plt.xticks(rotation=70, ha='right')  # Rotate labels by 45 degrees and align them to the right
    fig = plt.gcf()
    fig.set_size_inches(12, 8)  # Adjust the width and height as needed

    plt.show()


def rfr_classifier_wandb(data):
    wandb.init(project="project", entity="sb9342")

    # Log hyperparameters
    config = wandb.config
    config.n_estimators = 100
    config.max_depth = 12
    config.min_samples_split = 3
    config.min_samples_leaf = 3

    # Assuming X, y, and feature_labels are defined earlier
    X = data.drop(columns=['crmsd', 'name'])
    
    y = data['crmsd']
    feature_labels = X.columns  # Get the feature labels

    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    train_accuracy_history = []
    test_accuracy_history = []

    # Initialize the Random Forest Classifier
    for i in range(1, 100):  # Train for 100 iterations
        rfc = RandomForestClassifier(n_estimators=i, max_depth=12, min_samples_split=5, min_samples_leaf=5,)
        rfc.fit(X_train, y_train)

        # Calculate train and test accuracy
        train_pred = rfc.predict(X_train)
        test_pred = rfc.predict(X_test)
        train_accuracy = accuracy_score(y_train, train_pred)
        test_accuracy = accuracy_score(y_test, test_pred)

        # Append accuracy scores to history lists
        train_accuracy_history.append(train_accuracy)
        test_accuracy_history.append(test_accuracy)

        # Log metrics to W&B
        wandb.log({"Train Accuracy": train_accuracy, "Test Accuracy": test_accuracy})

    # Plot the change in test and train accuracy as the model runs
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, 100), train_accuracy_history, label='Train Accuracy', color='blue')
    plt.plot(range(1, 100), test_accuracy_history, label='Test Accuracy', color='red')
    plt.xlabel('Number of Trees')
    plt.ylabel('Accuracy')
    plt.title('Change in Test and Train Accuracy')
    plt.legend()
    plt.show()

    # Display feature importance
    feature_importance = pd.Series(rfc.feature_importances_, index=feature_labels)  # Use feature labels as index
    feature_importance.plot(kind='bar')
    plt.xlabel('Feature')
    plt.ylabel('Importance')
    plt.title('Feature Importance')
    plt.show()
    
def svm_classifier(data):
    X = data.drop(columns=['crmsd', 'name'])
    y = data['crmsd']

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    train_accuracy_history = []
    test_accuracy_history = []
    
    # Initialize the SVM Classifier
    for i in range(1, 50):  # Train for 50 iterations
        svm = SVC(kernel='linear', C=0.1*i)  # Linear kernel SVM with varying regularization parameter C
        svm.fit(X_train, y_train)

        # Calculate train and test accuracy
        train_pred = svm.predict(X_train)
        test_pred = svm.predict(X_test)
        train_accuracy = accuracy_score(y_train, train_pred)
        test_accuracy = accuracy_score(y_test, test_pred)

        # Append accuracy scores to history lists
        train_accuracy_history.append(train_accuracy)
        test_accuracy_history.append(test_accuracy)

    # Print the final train and test accuracy
    print("Final train accuracy:", train_accuracy)
    print("Final test accuracy:", test_accuracy)
    
    # Plot the change in test and train accuracy as the model runs
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, 50), train_accuracy_history, label='Train Accuracy', color='blue')
    plt.plot(range(1, 50), test_accuracy_history, label='Test Accuracy', color='red')
    plt.xlabel('Regularization Parameter (C)')
    plt.ylabel('Accuracy')
    plt.title('Change in Test and Train Accuracy')
    plt.legend()
    plt.show()

from sklearn.decomposition import PCA

def pca(data):
    print("helo")
    print(type(data))
    
    X = data.drop(columns=['crmsd', 'name', 'std_dev_1', 'std_dev2'])
    y = data['crmsd']
    feature_labels = X.columns  # Get the feature labels

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Initialize PCA and fit it to the training data
    pca = PCA(2)  # Choose number of components for visualization
    X_train_pca = pca.fit_transform(X_train)

    # Initialize the Random Forest Classifier
    train_accuracy_history = []
    test_accuracy_history = []
    for i in range(1, 100):  # Train for 100 iterations
        rfc = RandomForestClassifier(n_estimators=i, max_depth=12, min_samples_split=3, min_samples_leaf=3)
        rfc.fit(X_train_pca, y_train)

        # Calculate train and test accuracy
        train_pred = rfc.predict(X_train_pca)
        test_pred = rfc.predict(pca.transform(X_test))
        train_accuracy = accuracy_score(y_train, train_pred)
        test_accuracy = accuracy_score(y_test, test_pred)

        # Append accuracy scores to history lists
        train_accuracy_history.append(train_accuracy)
        test_accuracy_history.append(test_accuracy)

    # Plot the change in test and train accuracy as the model runs
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, 100), train_accuracy_history, label='Train Accuracy', color='blue')
    plt.plot(range(1, 100), test_accuracy_history, label='Test Accuracy', color='red')
    plt.xlabel('Number of Trees')
    plt.ylabel('Accuracy')
    plt.title('Change in Test and Train Accuracy')
    plt.legend()
    plt.show()

    # Display feature importance
    feature_importance = pd.Series(rfc.feature_importances_, index=range(pca.n_components_))  # Use PCA components as index
    feature_importance.plot(kind='bar')
    plt.xlabel('Principal Component')
    plt.ylabel('Importance')
    plt.title('Feature Importance (PCA)')
    plt.show()

import seaborn as sns
import shap

def rfr_feature(data):
    X = data.drop(columns=['crmsd', 'name'])
    y = data['crmsd']
    feature_labels = X.columns  # Get the feature labels

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
   
    train_accuracy_history = []
    test_accuracy_history = []
    # Initialize lists to store feature importance scores
    feature_importance_history = []
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'Avenir'
    # Set the color palette to green
    green_palette = ["#2ecc71"]  # Specify the shade of green you want
    sns.set_palette(green_palette)
# Set the font family to Times New Roman
    # Initialize the Random Forest Classifier
    for i in range(1, 100):  # Train for 100 iterations
        rfc = RandomForestClassifier(n_estimators=i, max_depth=12, min_samples_split=5, min_samples_leaf=5,)
        rfc.fit(X_train, y_train)
        train_pred = rfc.predict(X_train)
        test_pred = rfc.predict(X_test)
        train_accuracy = accuracy_score(y_train, train_pred)
        test_accuracy = accuracy_score(y_test, test_pred)

        # Append accuracy scores to history lists
        train_accuracy_history.append(train_accuracy)
        test_accuracy_history.append(test_accuracy)

        # Calculate feature importance
        #feature_importance = pd.Series(rfc.feature_importances_, index=feature_labels)
        #feature_importance_history.append(feature_importance)
    '''feature_importance.plot(kind='bar')
    plt.xlabel('Feature')
    plt.ylabel('Importance')
    plt.title('Feature Importance')
    plt.show() '''
    from pdpbox import pdp, info_plots

    #_________________pdp_______________________________________#
    def pred_func(rfc, X):
        predictions = rfc.predict_proba(X)

        if predictions.ndim == 1:
            predictions = predictions.reshape(-1, 1)
        return predictions
        
    print("train", train_accuracy)
    print("test", test_accuracy)
    # Calculate partial dependence for a specific feature (e.g., 'feature_name')
    ''' plot = info_plots.PredictPlot(
    df=X_test,
    feature='Predicted Hinge Count',
    feature_name='Predicted Hinge Count',
    pred_func= pred_func,
    model=rfc,
    model_features=feature_labels,
    grid_type='percentile',  # Adjust as needed
    num_grid_points=10       # Adjust the number of grid points as necessary
    )
    fig, axes, summary_df = plot.plot()
    print(fig)
    print(axes)
    print(summary_df)
    summary_df.to_csv('summary_data.csv', index=False)
    plt.show()
    '''
        #_________________pdp_______________________________________#

        #_________________pdp_______________________________________#

    
    # Calculate SHAP values
    explainer = shap.TreeExplainer(rfc)

    shap_values = explainer.shap_values(X_test)
    shap_sum = np.abs(shap_values).mean(axis=0)

    # Get the indices of the features sorted by importance
    importance_indices = np.argsort(shap_sum)[::-1]

    # Select the top N features
    top_n_indices = importance_indices[:10]  # Change 10 to the number of features you want to display
    top_features = X_test.columns[top_n_indices]

    # Plot SHAP values for the top features only
    shap.summary_plot(shap_values[:, top_n_indices], X_test.iloc[:, top_n_indices], feature_names=top_features)

    # Plot SHAP summary plot
    shap.summary_plot(shap_values, X_test, feature_names=['Hinge Count', 'Ratio Atoms'])
    plt.figure(figsize=(6, 16)) 
    plt.show()

    
# Generate the plot

    # Plot the change in test and train accuracy as the model runs
    
    ''' plt.figure(figsize=(10, 6))
    plt.plot(range(1, 100), train_accuracy_history, label='Train Accuracy', color='#1f77b4')

    # Plot smoothed test accuracy

    plt.plot(range(1, 100), test_accuracy_history, label='Test Accuracy', color='#ff7f0e')
    plt.xlabel('Number of Trees')
    plt.ylabel('Accuracy')
    plt.title('Change in Test and Train Accuracy')
    plt.legend()
    plt.show()'''

    # Plot the change in feature importance as the model runs
        
    
    # Plot the distribution of feature values for each label with enhanced aesthetics
     #sns.boxplot(x=y_train, y=X_train[feature], linewidth=2.5, fliersize=5)
        
        # Add labels and title with larger font size
    '''plt.xlabel('Label', fontsize=14)
        plt.ylabel(feature, fontsize=14)
        plt.title('Distribution of {} for Different Labels'.format(feature), fontsize=16)
        
        # Add grid lines for better readability
        plt.grid(True)
        # Show plot
        plt.show()'''

    # Analyze feature-label relationship for important features
    '''for feature in feature_importance.index:
        plt.figure(figsize=(10, 6))
        # Plot the distribution of feature values for each label
        sns.boxplot(x=y_train, y=X_train[feature])
        plt.xlabel('Label')
        plt.ylabel(feature)
        plt.title('Distribution of {} for Different Labels'.format(feature))
        plt.show()'''
    
if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Avenir'

# Set the style and color palette for the plots
    sns.set_style("whitegrid")
    sns.set_palette("pastel")
    data = preprocess('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/ml_analysis/all_features2.csv')
    rfr_feature(data)