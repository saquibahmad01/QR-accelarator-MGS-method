function [Q, R] = myQR(A)
%MYQR Summary of this function goes here
% The myQR function performs the following steps in to compute the QR
% decomposition of input matrix A based on Modified Gram Schmidt Algorithm using following steps:
% 1) Compute the euclidean norm of first column of A (a1)(This will be the
% diagonal entry of matrix R ) r11 =  ||a1||
% 2) Divide the first column of A (a1) by r11 to compute first column of Q
% (q1 = a1/r11)
% 3) The remaining row elements of R are calcuted using transpose of q1 and
% remaining A matrix (R12 = q1'A2)
% 4) Now update A using the following relation A' = A2 - q1*R12
% 5) Compute the QR factorization of the updated matrix A' 


%   Detailed explanation goes here
% The above function decomposes A in the following manner : 
% 1) The matrix R is computed in a row vector manner, and the matrix R is a
% upper traingular matrix which implies the elements below the main diagonal are all 0.  
% Further each diagonal entry of R is essentially the norm of the first
% column of A at a given iteration. The remaining elements of the row are
% the computed as a product of remaining A matrix and transpose of column q.
% 2) The matrix Q is computed in column vector manner where each column of
% q is essentially first column of A divided by the norm value at that
% stage
% 3) At the end of every stage the matrix A is updated and reduces by
% column in size. The updated matrix A which is computed by the above
% formula is now used for the QR decomposition 
% 4) This process continues until A is only one column wide. 

m = size(A,1);
n = size(A,2);

Q = zeros(m,n);
R = zeros(n,n);

q = zeros(m,1);

for i = 1:n 
    r = norm(A(:,1));
    q = A(:,1)/r;
    R_12 = q'*A(:,2:size(A,2));
    A = A(:,2:size(A,2)) - q*R_12;
    
    Q(:,i) = q;
    R(i,i) = r;
    R(i,i+1:n) = R_12;
end


end
