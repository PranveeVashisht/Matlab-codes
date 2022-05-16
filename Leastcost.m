format short
clear all
clc
%% input data
Cost = [2 10 4 5; 6 12 8 11; 3 9 5 7];
A = [12 25 20];  %row
B = [25 10 15 5];  %column

%% check balanced or unbalanced problem
if sum(A)==sum(B)
    fprintf('Given transportation problem is balanced\n');
else
    fprintf('Given transportation problem is unbalanced\n')
    if sum(A)<sum(B)
        Cost(end+1,:)= zeros(1,size(A,2));
        A(end+1)= sum(B)-sum(A);
    elseif sum(A)>sum(B)
        Cost(:,end+1) = zeros(1,size(A,2));
        B(end+1) = sum(A)- sum(B);
    end
end
ICost = Cost;  %save the cost copy
X = zeros(size(Cost));  %initialize allocation
[m,n] = size(Cost);  %Finding rows-column
BFS = m+n-1;  %Total BFS\

for  i = 1:size(Cost,1) %row loop
    for j = 1: size(Cost,2) %col loop
hh = min(Cost(:)); %finding minimum cost value
[rowind,colind] = find(hh==Cost);  %finding position of minimum cost cell

x11 = min(A(rowind), B(colind));  %assign allocations to each cost
[val,ind] = max(x11);  %finding max allocation
ii = rowind(ind);      %identify row posn
jj = colind(ind);      %identify col posn
y11 = min(A(ii), B(jj)); %find the value
X(ii,jj) = y11;  %assign ALLOCATION
A(ii) = A(ii)- y11; %reduce row value
B(jj) = B(jj)- y11; %reduce col value
Cost(ii,jj) = Inf;  %cell covered
    end
end
%% Print the initial BFS
fprintf('Initial BFS = \n');
IB = array2table(X);
disp(IB);
%% Check for degenerate and non-degenerate
TotalBFS = length(nonzeros(X));
if TotalBFS == BFS
    fprintf('Initial BFS is non-degenerate\n');
else
    fprintf('Initial BFS is degenerate\n');
end

%% Compute the initial transportation cost

InitialCost = sum(sum(ICost.*X));
fprintf('Initial BFS Cost = %d\n', InitialCost);









        

