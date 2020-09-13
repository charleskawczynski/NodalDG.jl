"""
function [F] = Filter1D(grid,Nc,s)
Purpose : Initialize 1D filter matrix of size N.
          Order of exponential filter is (even) s with cutoff at Nc;

# defaults based on suggestion, page 132
"""
function Filter1D(grid,Nc=0,s=16)

N = poly_order(grid)
filterdiag = ones(FT, N+1)
α = -log(eps(FT))

# Initialize filter function
for i=Nc:N
    filterdiag[i+1] = exp(-α*((i-Nc)/(N-Nc))^s);
end;

F = grid.V*diag(filterdiag)*grid.V⁻¹;
return F
end
