function weightMatrix = ccaAlgorithm(regulators,targets,p,n)

% regulators - matrix of data (methylation, copy number or expression) to
% represent potential regulators where rows are samples and columns are
% genes

% targets - matrix of data (always expression in our study) to represent
% potential targets - this must be the same size as the regulators matrix

% p - sampling subset size, suggested size 3-5

% n - number of iterations, suggested number depends on number of genes
% being tested. Check for convergence of the weight matrix

[a,b] = size(regulators)

weightMatrix = zeros(b,b);
regulatorSelected = zeros(b,1);
pairSelected = zeros(b,b);


for k = 1:n
 
    %select a set of genes of size p
    resample = randsample(b,p*2);
    subsample = resample(1:p);
    subsample2 = resample(p+1:end);
    
    
    targetSubsample = targets(:,subsample);
    regulatorSubsample = regulators(:, subsample2);
    
    [A,B,r,U,V] = canoncorr(targetSubsample,regulatorSubsample);
    
    C = abs(A);
    D = abs(B);
             
    Y = prctile(abs(B(:,1)),80);
    for i=1:p
        %regulators are i (row)
        %targets are j (cols)
        row = subsample2(i);
        %if regulator is greater than 80th percentile - score in relation
        %to the targets
        if abs(B(i,1)) > Y
           regulatorSelected(row,1) = regulatorSelected(row,1) + 1;
           for j=1:p
                col = subsample(j);              
                if row ~= col
                    pairSelected(row,col) = pairSelected(row,col) + 1;
                    weightMatrix(row,col) = weightMatrix(row,col) + ( (abs(A(j,1)) + abs(B(i,1)))*r(1))./((max(abs(A(:,1)))) + max(abs(B(:,1))));
                end
           end
        end
   end
    
end



for i =1:b
    for j=1:b
        if i ~= j
            weightMatrix(i,j) = ((weightMatrix(i,j)) / pairSelected(i,j));
        else
            weightMatrix(i,j) = .5;
        end
    end
end


end