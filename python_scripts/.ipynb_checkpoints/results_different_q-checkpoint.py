import os
import pandas as pd
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from scipy import stats
import matplotlib.pyplot as plt
import random

# --- Confidence Interval ---
def confidence_interval(accuracy, n, confidence=0.95):
    """
    Calculate the confidence interval for a given accuracy.

    :param accuracy: Accuracy value (float).
    :param n: Sample size (int).
    :param confidence: Confidence level (default 0.95).
    :return: Tuple (lower bound, upper bound).
    """
    std_error = np.sqrt(accuracy * (1 - accuracy) / n)
    z_score = stats.norm.ppf((1 + confidence) / 2)
    margin_of_error = z_score * std_error
    return (max(0, accuracy - margin_of_error), min(1, accuracy + margin_of_error))

# --- Load Data ---
def load_data(results_dir, network_name, simple_betas, complex_thetas):
    """
    Load the saved CSV files for simple and complex contagion simulations.

    :param results_dir: Directory where the results are saved.
    :param network_name: Name of the network.
    :param simple_betas: List of beta values for simple contagion.
    :param complex_thetas: List of theta values for complex contagion.
    :return: DataFrames for simple and complex contagion.
    """
    simple_data = []
    complex_data = []

    # Load simple contagion data
    for beta in simple_betas:
        file_path = os.path.join(results_dir, f'unweighted_simple_{network_name}_beta_{beta}_EPH.csv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            df['contagion_type'] = 'simple'
            df['beta'] = beta
            simple_data.append(df)

    # Load complex contagion data
    for theta in complex_thetas:
        for prob_complex in [-1,0.02, 0.03, 0.04]:
#         random.seed(theta)
#         prob_complex = random.choice([0.02, 0.03, 0.04])
#         prob_complex = random.choice([-1])
            file_path = os.path.join(results_dir, f'unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}_EPH-v2.csv')
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                df['contagion_type'] = 'complex'
                df['theta'] = theta
                df['prob_complex'] = prob_complex
                complex_data.append(df)
            else:
                print('DO NOT EXIST')
                print('don\'t exist')

    simple_df = pd.concat(simple_data, ignore_index=True) if simple_data else pd.DataFrame()
    complex_df = pd.concat(complex_data, ignore_index=True) if complex_data else pd.DataFrame()
    
    simple_df = simple_df.sample(len(simple_df),random_state=42)
    complex_df = complex_df.sample(len(complex_df),random_state=42)

    return simple_df, complex_df

# --- Decision Tree Classifier ---
def evaluate_classifier(simple_df, complex_df, thresholds, feature, random_state=42):
    """
    Train and evaluate a decision tree classifier.

    :param simple_df: DataFrame for simple contagion.
    :param complex_df: DataFrame for complex contagion.
    :param thresholds: List of theta values for complex contagion.
    :param feature: Feature to use ('EPH' or 'corr').
    :param random_state: Random state for reproducibility.
    :return: List of accuracies and confidence intervals for each threshold.
    """
    accuracies = []
    confidence_intervals = []

    for prob_complex in [-1,0.02, 0.03, 0.04]:
        # Filter complex contagion data for the current threshold
        complex_filtered = complex_df[complex_df['prob_complex'] == prob_complex]

        # Sample matching number of simple contagion data with random beta
        simple_sampled = simple_df.copy(deep=True)
        if len(simple_sampled) > len(complex_filtered):
            simple_sampled = simple_df.sample(n=len(complex_filtered), random_state=random_state)
        else:
            complex_filtered = complex_filtered.sample(n=len(simple_df), random_state=random_state)

        
        # Combine data and shuffle
        combined = pd.concat([simple_sampled, complex_filtered], ignore_index=True)
        combined = combined.sample(frac=1, random_state=random_state).reset_index(drop=True)

        combined = combined[~combined['EPH'].isna()]
        # Prepare features and labels
        if feature != 'merge':
            X = combined[[feature]]
        else:
            X = combined[["corr","EPH"]]
        
        y = combined['contagion_type'].map({'simple': 0, 'complex': 1})

        # Split into training and testing sets
        split_idx = int(0.7 * len(X))
        X_train, X_test = X[:split_idx], X[split_idx:]
        y_train, y_test = y[:split_idx], y[split_idx:]

        # Train decision tree
        clf = DecisionTreeClassifier(random_state=random_state)
        clf.fit(X_train, y_train)

        # Evaluate accuracy
        y_pred = clf.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        accuracies.append(accuracy)

        # Calculate confidence interval
        ci = confidence_interval(accuracy, len(y_test))
        confidence_intervals.append(ci)

    return accuracies, confidence_intervals

# --- Plot Accuracy ---
def plot_accuracy(thresholds, accuracies, confidence_intervals, feature, output_file):
    """
    Plot accuracy vs. threshold for a given feature.

    :param thresholds: List of thresholds.
    :param accuracies: List of accuracies for each threshold.
    :param confidence_intervals: List of confidence intervals for each threshold.
    :param feature: Feature used ('EPH' or 'corr').
    :param output_file: Path to save the figure.
    """
    lower_bounds = [ci[0] for ci in confidence_intervals]
    upper_bounds = [ci[1] for ci in confidence_intervals]

    plt.figure(figsize=(8, 6))
    plt.plot(thresholds, accuracies, marker='o', label=feature, linestyle='-')
    plt.fill_between(thresholds, lower_bounds, upper_bounds, alpha=0.2)
    plt.xlabel('q', fontsize=12)
    plt.ylabel('Accuracy', fontsize=12)
    plt.title(f'Accuracy vs. Threshold using {feature}', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()
    print(f"Figure saved to {output_file}")
    
def plot_accuracy_together(thresholds,  accuracies_eph, accuracies_corr, accuracies_merge
                      , confidence_intervals_eph, confidence_intervals_corr, confidence_intervals_merge, output_file):
    
    lower_bounds_eph = [ci[0] for ci in confidence_intervals_eph]
    upper_bounds_eph = [ci[1] for ci in confidence_intervals_eph]
    
    lower_bounds_corr = [ci[0] for ci in confidence_intervals_corr]
    upper_bounds_corr = [ci[1] for ci in confidence_intervals_corr]
    
    lower_bounds_merge = [ci[0] for ci in confidence_intervals_merge]
    upper_bounds_merge = [ci[1] for ci in confidence_intervals_merge]
    
    # Create the plot
    plt.figure(figsize=(8, 6))

    plt.plot([0,0.02, 0.03, 0.04], accuracies_eph, marker='o', label='EPH')
    plt.fill_between([0,0.02, 0.03, 0.04], lower_bounds_eph, upper_bounds_eph, alpha=0.2)

    plt.plot([0,0.02, 0.03, 0.04], accuracies_corr, marker='s', label='corr(order,deg)')
    plt.fill_between([0,0.02, 0.03, 0.04], lower_bounds_corr, upper_bounds_corr, alpha=0.2)
    
#     plt.plot(thresholds, accuracies_merge, marker='p', label='corr & EPH')
#     plt.fill_between(thresholds, lower_bounds_merge, upper_bounds_merge, alpha=0.2)
    
    # Set x-axis to show only integer values
    plt.xticks([0,0.02, 0.03, 0.04])


    plt.xlabel('q',fontsize=15)
    plt.ylabel('Accuracy',fontsize=15)
    plt.legend(loc='lower right',fontsize=15)    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_file, format='pdf', dpi=1200, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory
    
    
    

# --- Main Function ---
def main():
    # Parameters
    ['email_eu_modified', 'conf','hospital', 'school', 'work']
    for network_name in ['email_eu_modified', 'conf', 'school']:
        results_dir = '../result_85stop_10init_nodes-v2/unweighted/'
        simple_betas = [0.02, 0.03, 0.04, 0.05]
        complex_thetas = [2, 3, 4, 5, 6]
        
        #both together
        simple_df, complex_df = load_data(results_dir, network_name, simple_betas, complex_thetas)
        
        accuracies_eph, confidence_intervals_eph = evaluate_classifier(simple_df, complex_df, complex_thetas, 'EPH')
        accuracies_corr, confidence_intervals_corr = evaluate_classifier(simple_df, complex_df, complex_thetas, 'corr')
        accuracies_merge, confidence_intervals_merge = evaluate_classifier(simple_df, complex_df, complex_thetas, 'merge')
        
        output_file = f'../result_85stop_10init_nodes-v2/accuracy_vs_q_{network_name}-v2.pdf'

        plot_accuracy_together(complex_thetas, accuracies_eph, accuracies_corr, accuracies_merge
                      , confidence_intervals_eph, confidence_intervals_corr, confidence_intervals_merge, output_file)
        

#         # Load data
#         simple_df, complex_df = load_data(results_dir, network_name, simple_betas, complex_thetas)

#         # Evaluate classifiers for both features
#         for feature in ['EPH', 'corr']:
#             accuracies, confidence_intervals = evaluate_classifier(simple_df, complex_df, complex_thetas, feature)

#             # Plot and save figure
#             output_file = f'../results/accuracy_vs_threshold_{feature}_{network_name}.png'
#             plot_accuracy(complex_thetas, accuracies, confidence_intervals, feature, output_file)

if __name__ == "__main__":
    main()
