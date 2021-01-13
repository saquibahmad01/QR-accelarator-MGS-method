function [MSE] = getMSE(Q, R, A)
%GETMSE Summary of this function goes here
%   Detailed explanation goes here
for i = 1:size(A,3)
   Qerror = Q(:,:,i)'*Q(:,:,i) - eye(size(A,2),size(A,2));
   Aerror = Q(:,:,i)*R(:,:,i) - A(:,:,i);
   MSEQ(i) = norm(Qerror(:)).^2/length(Qerror(:));
   MSEA(i) = norm(Aerror(:)).^2/length(Qerror(:));
end

MSE = mean([MSEQ, MSEA]);

end

