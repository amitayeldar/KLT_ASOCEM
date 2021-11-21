function g = NeumannBoundCondMod(f,h)
% Neumann boundary condition for derivative of step size h

%%
switch nargin
    case 1
        h = 1;
    case 2
    otherwise
end
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([1+h nrow-h],[1+h ncol-h]);  
g([1 nrow],2:end-1) = g([1+h nrow-h],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[1+h ncol-h]);  
end