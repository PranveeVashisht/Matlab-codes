% multi-objective LPP
clc
clear

A = [2 4 1 0 1 0; 3 5 4 -1 0 1];
b = [8; 15];
C = [3 2 4 0 0; 1 5 3 0 0]
C1 = [0 0 0 0 0 1];     % phase-I objective function

n_art = 1;
[m, n] = size(A);

C = sum(C)/m;
A = [A b];

% Phase-I
fprintf('* PHASE-1 * \n')

cost = zeros(1, n + 1);
cost(1:n) = C1;

bv = (n - m + 1):(n);   % basic variables
zjcj = cost(bv)*A - cost;
zcj = [zjcj; A];

smpTb=array2table(zcj);
smpTb.Properties.VariableNames(1:n+1) = ["x1", "x2", "x3", "s2", "s1", "a2", "sol"];
disp(smpTb)

RUN=true;
while RUN

    if any(zjcj(1:n)>0) 
         fprintf(' the current BFS is not optimal \n')
         zc=zjcj(1:n);

         [Enter_val, pvt_col]= max(zc);

         if all(A(:,pvt_col)<0)
            error('LPP is Unbounded' ,pvt_col);
         else
            sol=A(:,end);
            column=A(:,pvt_col);
             for i=1:size(A)
                if column(i)>0
                   ratio(i)= sol(i)./column(i);
                else
                   ratio(i)=inf;
               end
           end
           [leaving_val, pvt_row]=min(ratio);
        end
        bv(pvt_row)=pvt_col;
        pvt_key=A(pvt_row, pvt_col);
        A(pvt_row,:)=A(pvt_row,:)./pvt_key;
        for i=1:size(A,1)
             if i~=pvt_row
                A(i,:)=A(i,:)-A(i, pvt_col).*A(pvt_row,:);
             end
        end

        zjcj = zjcj - zjcj(pvt_col).*A(pvt_row, :);
        zcj=[zjcj;A];
        smpTb=array2table(zcj);
        smpTb.Properties.VariableNames(1:n+1) = ["x1", "x2", "x3", "s2", "s1", "a2", "sol"];
        disp(smpTb)
    else
        RUN=false;
        disp("Current BFS is optimal.")
         
    end
end

fprintf('* PHASE-2 * \n')
cost = zeros(1, n - n_art + 1);
cost(1: n - n_art) = C;

A(:, n - n_art + 1: n) = [];
zjcj = cost(bv)*A - cost;
zcj = [zjcj; A];

n = n - n_art;

smpTb=array2table(zcj)
smpTb.Properties.VariableNames(1:n+1) = ["x1", "x2", "x3", "s2", "s1", "sol"];
disp(smpTb)

RUN=true;
while RUN
    if any(zjcj(1:end-1)<0)
        fprintf(' the current BFS is not optimal \n')
        zc=zjcj(1:end-1);
        [Enter_val, pvt_col]= min(zc);
        if all(A(:,pvt_col)<=0)
            error('LPP is Unbounded all enteries are <=0 in column %d',pvt_col);
        else
            sol=A(:,end);
            column=A(:,pvt_col);
            for i=1:size(A,1)
                if column(i)>0
                    ratio(i)= sol(i)./column(i);
                else
                    ratio(i)=inf;
                end
            end
            [leaving_val, pvt_row]=min(ratio);
        end
        bv(pvt_row)=pvt_col;
        pvt_key=A(pvt_row, pvt_col);
        A(pvt_row,:)=A(pvt_row,:)./pvt_key;
        for i=1:size(A,1)
             if i~=pvt_row
                A(i,:)=A(i,:)-A(i, pvt_col).*A(pvt_row,:);
             end
        end
        zjcj = zjcj - zjcj(pvt_col).*A(pvt_row, :);
        zcj=[zjcj;A];
        smpTb=array2table(zcj);
        smpTb.Properties.VariableNames(1:n+1) = ["x1", "x2", "x3", "s2", "s1", "sol"];
        disp(smpTb)
    else
        RUN=false;
        fprintf('The current BFS is optimal \n');
        z=input('Enter 0 for minimization and 1 for max \n');
        if z==0
            Obj_value=-zjcj(end);
        else
            Obj_value=zjcj(end);
        end
        fprintf('The final optimal value is %f\n',Obj_value);
    end
end