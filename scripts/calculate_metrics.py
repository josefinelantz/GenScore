def calculate_metrics(y_true, y_pred, y_scores):
    metrics = {
        "Precision": precision_score(y_true, y_pred, average="binary"),
        "Recall": recall_score(y_true, y_pred, average="binary"),
        "F1-Score": f1_score(y_true, y_pred, average="binary"),
        "AUC": roc_auc_score(y_true, y_scores)
    }
    return metrics

def plot_roc_curve(y_true, y_scores):
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    plt.plot(fpr, tpr, label="ROC Curve")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend()
    plt.show()