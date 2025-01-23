function C=FlucDis(J,D)
% Fluctuation-Dissipation Theorem
% JC+CJ^T+D=0
% J is Jacobian Matrix
% C is Correlation Matrix
% D is Diffusion Matrix

% Sylvester Equation AX+XB=C
% In our case, A=J, B=J^T, C=-D, and X=C.
% C=sylvester(J,J',-D) works numerically.
% However, we can use the Kronecker product notation and the vectorization operator to
% try and solve analytically.

% Examples:
% https://math.stackexchange.com/questions/3144254/solve-symbolic-sylvester-like-equation-in-matlab-or-maple
% https://math.stackexchange.com/questions/353329/solution-of-a-sylvester-equation

vecD=D(:); % A vector representing the diffusion matrix.
vecD=-vecD; % On the other side of the equal sign.
I=eye(size(J,1)); % Create and identity matrix equivalent to the size of the Jacobian.
K=kron(I,J)+kron(J,I);
% K*vecC=vecD implies vecC=inv(K)*vecD=adjoint(K)/det(K)*vecD
% vecC=inv(K)*vecD;
vecC=adjoint(K)/det(K)*vecD;
C=reshape(vecC,[size(J,1),size(J,1)]);
end