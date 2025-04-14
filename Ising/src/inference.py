import numpy as np
import torch
from torch import optim
from torch.optim.lr_scheduler import ExponentialLR
from tqdm import tqdm

def gamma(J0, J1):
    """
    Computes the normalized Frobenius norm of the difference between two coupling matrices.

    Parameters:
        J0 (ndarray): The true (hidden) coupling matrix.
        J1 (ndarray): The inferred coupling matrix.

    Returns:
        float: The normalized Frobenius norm of the difference.
    """
    value = np.sqrt(np.sum((J0 - J1) ** 2) / np.sum(J0 ** 2))
    return value


def evaluate_MF_inference(spins_matrix, beta):
    """
    Performs mean-field inference to estimate the coupling matrix using the provided spin data.

    Parameters:
        spins_matrix (ndarray): Data matrix with shape (samples, variables), where each row is a spin configuration.
        beta (float): Inverse temperature parameter.

    Returns:
        ndarray: The inferred coupling matrix. If the covariance matrix is singular, returns a matrix filled with NaN.
    """
    # Compute the covariance matrix of the spin data
    cov_matrix = np.cov(spins_matrix.T)
    try:
        # Attempt to compute the inverse of the covariance matrix
        cov_inv = np.linalg.inv(cov_matrix)
        # Set the diagonal elements to zero (no self-coupling)
        np.fill_diagonal(cov_inv, 0)
        # Return the inferred coupling matrix scaled by -1/beta
        return -cov_inv / beta
    except np.linalg.LinAlgError:
        # If the covariance matrix is singular, return a matrix filled with NaN
        error_matrix = np.zeros_like(cov_matrix)
        error_matrix[:] = np.nan
        return error_matrix
    
def evaluate_PL_inference(data, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.0):
    """
    Performs pseudo-likelihood inference to estimate the coupling matrix
    using the provided data. All auxiliary functions are nested inside.

    Parameters:
        data (array-like): Data matrix with shape (samples, variables).
        learning_rate (float): Initial learning rate for optimization.
        epochs (int): Number of training epochs.
        decay (float): Exponential decay rate of the learning rate.
        l1_amplitude (float): Amplitude for the L1 regularization term.

    Returns:
        final_couplings (ndarray): The inferred coupling matrix with zeros on the diagonal.
    """
    # Convert input data to a PyTorch tensor
    data = torch.Tensor(data)
    # Assume the number of variables is given by the number of columns
    N = data.shape[1]

    def singleRow_loss(couplings, batch, row_index):
        """
        Computes the negative log-likelihood loss for a single variable (row)
        given the current coupling matrix.
        """
        # Compute the effective field for the current variable
        effective_field = torch.mm(couplings[row_index].view(1, -1), batch.t()).view(-1)
        # Remove self-coupling contribution
        effective_field = effective_field - couplings[row_index][row_index] * batch[:, row_index]
        # Compute the exponent part of the loss (with a factor -2*s_i*h_i)
        exponent = -2 * batch[:, row_index] * effective_field
        # Compute the logistic function to get probabilities
        loss_vector = 1 / (1 + torch.exp(exponent))
        # Calculate the negative log-likelihood for the current row
        loss = -torch.mean(torch.log(loss_vector))
        return loss

    def loss_function(couplings, batch, l1_amplitude):
        """
        Aggregates the losses computed for each variable and adds an L1 penalty
        for regularization.
        """
        total_loss = 0
        for row_index in range(N):
            total_loss += singleRow_loss(couplings, batch, row_index)
        # Average the total loss and add L1 regularization term to encourage sparsity
        return total_loss / N + l1_amplitude * torch.mean(torch.abs(couplings))

    # Initialize couplings with random values.
    # The variance is scaled with 1/N and diagonal elements are later set to zero.
    couplings = torch.normal(0, 1 / N, size=(N, N))
    couplings = couplings - torch.diag(couplings.diag())

    # Create an optimized couplings tensor with gradient tracking enabled
    optimized_couplings = couplings.clone().requires_grad_(True)

    # Setup the optimizer (SGD) and learning rate scheduler
    optimizer = optim.SGD([optimized_couplings], lr=learning_rate)
    scheduler = ExponentialLR(optimizer, gamma=decay)

    # Training loop over the specified number of epochs
    for epoch in tqdm(range(epochs), desc="Training Progress", unit="epoch"):
        # Compute the loss using the nested loss_function
        loss = loss_function(optimized_couplings, data, l1_amplitude)
        optimizer.zero_grad()  # Zero out any previous gradients
        loss.backward()        # Compute gradients
        optimizer.step()       # Update the coupling parameters
        scheduler.step()       # Decay the learning rate

    # Detach the final optimized couplings and convert to a NumPy array;
    # ensure that the diagonal elements remain zero.
    final_couplings = optimized_couplings.detach().numpy()
    np.fill_diagonal(final_couplings, 0)

    final_couplings /= beta  # Scale the couplings by the inverse temperature
    return final_couplings
    
