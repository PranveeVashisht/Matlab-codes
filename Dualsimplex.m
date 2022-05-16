%simplex method
% Max. z = x1 + 2x2, subject to
% − x1 + x2 ≤ 1,
% x1 + x2 ≤ 2,
% x1, x2 ≥ 0.
clc;
clear;
C = [-3 -5 0 0 0];
A = [-1 -3 1 0 -3;-1 -1 0 1 -2];
bv = [3 4];

[m,n] = size(A);

ZjCj= C(bv)*A - C;
ZCj= [ZjCj; A];
%bi = ZCj[]
Simplextable= array2table(ZCj);
Simplextable.Properties.VariableNames(1:size(Simplextable,2))= {'x1','x2','s1','s2','sol' };
fprintf('\n');
disp(Simplextable);

% matrix for storing values of x1, x2,...xn in final solution
var= zeros(1:n);

flag= 1;
while flag
    b=A(:,end);
    if any(b < 0)
        fprintf("The current BFS is not feasib;e\n");
        ZC=ZjCj(1:end-1);
        [val1 , pvt_row]= min(b);
        if all(A(pvt_row,:) >= 0)
            fprintf("LPP has unbounded solution");
        else
            ZC=ZjCj(1:end-1);
            for i=1:size(A,2)-1
                if(A(pvt_row,i)<0)
                    ratio(i)=abs(ZC(i)/A(pvt_row,i));
                else
                    ratio(i)=inf;
                end
                [val2 , pvt_col]= min(ratio);
            end
           

            bv(pvt_row)= pvt_col ;  %added 1 because A matrix begins from 2nd row in ZCj
            pivot= A(pvt_row, pvt_col);
            A(pvt_row,:) = A(pvt_row,:) ./ pivot;
           
            for i= 1: m
                if(i~= pvt_row)
                    A(i,:)= A(i,:) - (A(i,pvt_col)*A(pvt_row,:));
                end
            end

            ZjCj= ZjCj - (ZjCj(pvt_col).*A(pvt_row,:));  
            ZCj= [ZjCj; A];
            table=array2table(ZCj);
            table.Properties.VariableNames(1:size(Simplextable,2))= {'x1','x2','s1','s2','sol'};

        end
    else
        fprintf("The current BFS is optimal");
        disp(table)
        flag= 0;
    end
end