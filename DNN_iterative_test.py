import pandas as pd
import numpy as np
import torch
from torch import nn, optim
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from sklearn.metrics import mean_squared_error, mean_absolute_error

# Check if a GPU is available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Load data from specified path
data = pd.read_csv('/home/zhuy27/kv_cache/test/optimization/ML_data.csv', header=None)

# # Normalize each column
# scaler = MinMaxScaler()  # Use StandardScaler() for standardization
# data = scaler.fit_transform(data)

# Step 1: Prepare input-output pairs for training and testing
n_trajectories = len(data) // 83
train_input, train_output = [], []
test_input, test_target = [], []

# Split based on trajectories
train_trajectories = int(0.7 * n_trajectories)  # 70% for training
for i in range(n_trajectories):
    trajectory = data.iloc[i * 83:(i + 1) * 83].values

    # # Remove the fifth column from the trajectory data
    # trajectory = np.delete(trajectory, 4, axis=1)

    if i < train_trajectories:
        # Training: use all timesteps in the trajectory except the last
        train_input.extend(trajectory[:-1, :])    # All but the last timestep as input
        train_output.extend(trajectory[1:, :])    # All but the first timestep as target
    else:
        # Testing: only use the first timestep as input and the last timestep as target
        test_input.append(trajectory[0, :])       # First timestep as input
        test_target.append(trajectory[-1, :])     # Last timestep as target


# Convert training data to PyTorch tensors
train_input_tensor = torch.tensor(np.array(train_input), dtype=torch.float32).to(device)
train_output_tensor = torch.tensor(np.array(train_output), dtype=torch.float32).to(device)

# Convert testing data to PyTorch tensors
test_input_tensor = torch.tensor(np.array(test_input), dtype=torch.float32).to(device)
test_target_tensor = torch.tensor(np.array(test_target), dtype=torch.float32).to(device)

# Create DataLoader for training data
train_dataset = TensorDataset(train_input_tensor, train_output_tensor)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

# define the DNN layer parameters as arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--num_layers", type=int, default=2, help="Number of hidden layers in the DNN")
parser.add_argument("--hidden_size", type=int, default=128, help="Hidden layer size in the DNN")
args = parser.parse_args()

# Define the DNN Model for Single-Step Prediction
class SimpleDNN(nn.Module):
    def __init__(self, input_size=12, hidden_size=args.hidden_size, output_size=12, num_layers=args.num_layers):
        super(SimpleDNN, self).__init__()
        layers = []
        for i in range(num_layers):
            layers.append(nn.Linear(input_size if i == 0 else hidden_size, hidden_size))
            layers.append(nn.ReLU())
        layers.append(nn.Linear(hidden_size, output_size))
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

# Instantiate the model
model = SimpleDNN().to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Set update frequency for tqdm
update_frequency = 1000

# Training Loop with tqdm
n_epochs = 40
for epoch in range(n_epochs):
    model.train()
    train_loss = 0.0
    with tqdm(total=len(train_loader), desc=f"Epoch {epoch+1}/{n_epochs}", unit="batch") as tepoch:
        for batch_idx, (inputs, targets) in enumerate(train_loader):
            inputs, targets = inputs.to(device), targets.to(device)

            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            # Update tqdm only every `update_frequency` batches
            if (batch_idx + 1) % update_frequency == 0:
                tepoch.set_postfix(loss=loss.item())
                tepoch.update(update_frequency)
    
    train_loss /= len(train_loader)
    print(f"Epoch {epoch + 1}/{n_epochs}, Training Loss: {train_loss:.4f}")

# Save the model after training
model_path = f"/home/zhuy27/kv_cache/test/optimization/model_checkpoint/trained_dnn_model_{args.num_layers}_{args.hidden_size}.pth"
torch.save(model.state_dict(), model_path)
print(f"Model saved to {model_path}")

# Testing Phase with Iterative Prediction and Additional Metrics
model.eval()
all_predictions = []
all_targets = []
all_trace = []
with torch.no_grad():
    # Iterate over each test trajectory individually
    for i in range(len(test_input_tensor)):
        # Initialize current_input with the first timestep in the sequence
        current_input = test_input_tensor[i].unsqueeze(0)  # Shape: (1, features)
        all_trace.append(current_input.cpu().numpy())

        # Iterate 82 times to predict each timestep sequentially
        for _ in range(82):
            current_output = model(current_input)  # Predict next timestep
            current_input = current_output  # Set predicted output as the next input
            all_trace.append(current_input.cpu().numpy())
        
        # `current_output` is now the model's prediction for the last timestep
        current_output = current_output.view(-1)  # Reshape to match the target's shape
        all_predictions.append(current_output.cpu().numpy())
        all_targets.append(test_target_tensor[i].cpu().numpy())  # True last timestep target

# Calculate Metrics
all_predictions = np.concatenate(all_predictions, axis=0)
all_targets = np.concatenate(all_targets, axis=0)
all_trace = np.concatenate(all_trace, axis=0)

# Assuming `all_predictions` and `all_targets` are already populated as lists of 1D arrays (predictions per sample)
all_predictions = np.array(all_predictions).reshape(-1, 12)  # Reshape to have 11 columns
all_targets = np.array(all_targets).reshape(-1, 12)  # Reshape to have 11 columns
all_trace = np.array(all_trace).reshape(-1, 12)  # Reshape to have 11 columns

# Calculate Metrics
mse = mean_squared_error(all_targets, all_predictions)
mae = mean_absolute_error(all_targets, all_predictions)
rmse = np.sqrt(mse)

print(f"Test MSE: {mse:.4f}")
print(f"Test MAE: {mae:.4f}")
print(f"Test RMSE: {rmse:.4f}")

# Save all_targets and all_predictions to separate CSV files
targets_df = pd.DataFrame(all_targets, columns=[f"target_{i}" for i in range(all_targets.shape[1])])
predictions_df = pd.DataFrame(all_predictions, columns=[f"prediction_{i}" for i in range(all_predictions.shape[1])])
trace_df = pd.DataFrame(all_trace, columns=[f"trace_{i}" for i in range(all_trace.shape[1])])

# Save to CSV files
targets_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_targets_DNN_{args.num_layers}_{args.hidden_size}.csv", index=False)
predictions_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_predictions_DNN__{args.num_layers}_{args.hidden_size}.csv", index=False)
trace_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_trace_DNN__{args.num_layers}_{args.hidden_size}.csv", index=False)

print("Targets saved to test_targets.csv")
print("Predictions saved to test_predictions.csv")