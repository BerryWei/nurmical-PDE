import torch
import scipy



def lap2d(dx, dy, Nx, Ny):
    Ax = (1 / (dx * dx)) * (torch.diag(-2 * torch.ones(Nx)) + torch.diag(torch.ones(Nx - 1), diagonal=-1) +
                             torch.diag(torch.ones(Nx - 1), diagonal=1))
    Ay = (1 / (dy * dy)) * (torch.diag(-2 * torch.ones(Ny)) + torch.diag(torch.ones(Ny - 1), diagonal=-1) +
                             torch.diag(torch.ones(Ny - 1), diagonal=1))
    
    L = torch.kron(torch.eye(Ny), Ax) + torch.kron(Ay, torch.eye(Nx))
    
    return L

def torch_to_scipy_sparse(tensor: torch.tensor) -> scipy.sparse:


    # Get zero idx, value
    nonzero_indices = torch.nonzero(tensor, as_tuple=True)
    nonzero_values = tensor[nonzero_indices].cpu().numpy()

    sparse_matrix = scipy.csr_matrix((nonzero_values, nonzero_indices), shape=tensor.shape)

    return sparse_matrix