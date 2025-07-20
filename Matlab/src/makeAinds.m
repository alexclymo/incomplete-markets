function [rowInd_A_Xup,colInd_A_Xup,rowInd_A_Xdown,colInd_A_Xdown] = makeAinds(Na,Nz)

% A: upper and lower diagonals for asset drift. These are saved as a long
% column vector corresponding to the indices for each flattened x,z
% row. Start constructing as matrix then flatten
rowInd_A_Xup = zeros(Na-1,Nz);
colInd_A_Xup = zeros(Na-1,Nz);
rowInd_A_Xdown = zeros(Na-1,Nz);
colInd_A_Xdown = zeros(Na-1,Nz);
for i = 1:Nz %loop over current z
    % upper diagonal
    rowInd_A_Xup(:,i) = (1:Na-1)' + Na*(i-1);
    colInd_A_Xup(:,i) = (2:Na)' + Na*(i-1);
    % lower diagonal
    rowInd_A_Xdown(:,i) = (2:Na)' + Na*(i-1);
    colInd_A_Xdown(:,i) = (1:Na-1)' + Na*(i-1);
end
rowInd_A_Xup = rowInd_A_Xup(:);
colInd_A_Xup = colInd_A_Xup(:);
rowInd_A_Xdown = rowInd_A_Xdown(:);
colInd_A_Xdown = colInd_A_Xdown(:);