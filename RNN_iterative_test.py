import pandas as pd
import numpy as np
import torch
from torch import nn, optim
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from sklearn.metrics import mean_squared_error, mean_absolute_error


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


data = pd.read_csv('/home/zhuy27/kv_cache/test/optimization/ML_data.csv', header=None)


sequence_length = 25 
n_features = data.shape[1]  


n_trajectories = len(data) // 83
train_sequences, train_targets = [], []
test_sequences, test_targets = [], []

train_trajectories = int(0.7 * n_trajectories) 

for i in range(n_trajectories):
    trajectory = data.iloc[i * 83:(i + 1) * 83].values

    if i < train_trajectories:
        
        for j in range(len(trajectory) - sequence_length):
            seq_input = trajectory[j:j + sequence_length, :]  
            seq_output = trajectory[j + sequence_length, :]   
            train_sequences.append(seq_input)
            train_targets.append(seq_output)
    else:
        
        test_sequences.append(trajectory[:sequence_length, :])
        test_targets.append(trajectory[-1, :])


train_sequences_tensor = torch.tensor(np.array(train_sequences), dtype=torch.float32).to(device)  # (num_samples, seq_len, n_features)
train_targets_tensor = torch.tensor(np.array(train_targets), dtype=torch.float32).to(device)      # (num_samples, n_features)

test_sequences_tensor = torch.tensor(np.array(test_sequences), dtype=torch.float32).to(device)    # (num_test_samples, seq_len, n_features)
test_targets_tensor = torch.tensor(np.array(test_targets), dtype=torch.float32).to(device)        # (num_test_samples, n_features)


train_dataset = TensorDataset(train_sequences_tensor, train_targets_tensor)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)


class SimpleRNN(nn.Module):
    def __init__(self, input_size=12, hidden_size=128, output_size=12, num_layers=2):
        super(SimpleRNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        # LSTM
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True)

        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(device)  # (num_layers, batch, hidden_size)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(device)  # (num_layers, batch, hidden_size)
        

        out, _ = self.lstm(x, (h0, c0))  # out: (batch, seq_len, hidden_size)
        

        out = out[:, -1, :]  # (batch, hidden_size)
        

        out = self.fc(out)  # (batch, output_size)
        return out


model = SimpleRNN(input_size=n_features, hidden_size=128, output_size=n_features, num_layers=2).to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)


update_frequency = 1000


n_epochs = 40
for epoch in range(n_epochs):
    model.train()
    train_loss = 0.0
    with tqdm(total=len(train_loader), desc=f"Epoch {epoch+1}/{n_epochs}", unit="batch") as tepoch:
        for batch_idx, (inputs, targets) in enumerate(train_loader):
            # inputs: (batch, seq_len, n_features)
            # targets: (batch, n_features)
            inputs, targets = inputs.to(device), targets.to(device)

            optimizer.zero_grad()
            outputs = model(inputs)  # outputs: (batch, n_features)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
     
            if (batch_idx + 1) % update_frequency == 0:
                tepoch.set_postfix(loss=loss.item())
                tepoch.update(update_frequency)
            else:
                tepoch.update(1)
    
    train_loss /= len(train_loader)
    print(f"Epoch {epoch + 1}/{n_epochs}, Training Loss: {train_loss:.4f}")


model_path = f"/home/zhuy27/kv_cache/test/optimization/model_checkpoint/trained_rnn_model_{sequence_length}.pth"
torch.save(model.state_dict(), model_path)
print(f"Model saved to {model_path}")


model.load_state_dict(torch.load(model_path))



model.eval()
all_predictions = []
all_targets = []
all_trace = []

with torch.no_grad():
    for i in range(len(test_sequences_tensor)):
        
        current_input = test_sequences_tensor[i].unsqueeze(0)  # (1, seq_len, n_features)
        trace = current_input.cpu().numpy()

       
        for _ in range(83 - sequence_length):
            
            current_output = model(current_input)  # (1, n_features)
            
            
            current_input = torch.cat((current_input[:, 1:, :], current_output.unsqueeze(1)), dim=1)  # (1, seq_len, n_features)
            trace = np.concatenate((trace, current_output.cpu().numpy()[:, np.newaxis, :]), axis=1)
        
        all_predictions.append(current_output.cpu().numpy().flatten())
        all_targets.append(test_targets_tensor[i].cpu().numpy())
        all_trace.append(trace)


all_predictions = np.array(all_predictions)  # (num_test_samples * (83 - seq_len), n_features)
all_targets = np.array(all_targets)          # (num_test_samples, n_features)
all_trace = np.array(all_trace)              # (num_test_samples, 83, n_features)


final_predictions = all_predictions
final_targets = all_targets


mse = mean_squared_error(final_targets, final_predictions)
mae = mean_absolute_error(final_targets, final_predictions)
rmse = np.sqrt(mse)

print(f"Test MSE: {mse:.4f}")
print(f"Test MAE: {mae:.4f}")
print(f"Test RMSE: {rmse:.4f}")


targets_df = pd.DataFrame(final_targets, columns=[f"target_{i}" for i in range(final_targets.shape[1])])
predictions_df = pd.DataFrame(final_predictions, columns=[f"prediction_{i}" for i in range(final_predictions.shape[1])])
trace_df = pd.DataFrame(all_trace.reshape(-1, all_trace.shape[2]), columns=[f"trace_{i}" for i in range(all_trace.shape[2])])


targets_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_targets_RNN_{sequence_length}.csv", index=False)
predictions_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_predictions_RNN_{sequence_length}.csv", index=False)
trace_df.to_csv(f"/home/zhuy27/kv_cache/test/optimization/result/test_trace_RNN_{sequence_length}.csv", index=False)

print("Targets saved to test_targets.csv")
print("Predictions saved to test_predictions.csv")
print("Trace saved to test_trace.csv")
