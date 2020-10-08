function [solution,errorConvergenceRate,numOfGSIterations] = GSFunction1(globalA, globalB, height, numOfNodes, inputUVector, meshNum)
% For "Test GS1": add "errorConvergenceRate,numOfGSIterations" as
%                 more function outputs.
% For "With Multigrid": Runs for 1 iteration.


    U_i               = inputUVector;
    
    % Test GS
    normU_iVector     = zeros(1,1);
    errors2           = zeros(1,1);
    norm_0            = getHrCurlNormforProblem3(inputUVector,globalA);
    checkForTolerance = 1;
    numOfGSIterations = 0;
    
    
    mystr             = ['ZMatrix/ZMatrix' num2str(meshNum) '.mat'];
    load(mystr);
    

    
    while checkForTolerance >= (10^(-15))             % For "Test GS".
        numOfGSIterations    = numOfGSIterations +1;  % For "Test GS".
        
        for i = 1:numOfNodes
            
            
            Z_i                                = ZMatrix(i,:);
            Z_i                                = Z_i(Z_i > 0);
            diff                               = globalB - globalA * U_i;
            diff                               = diff(Z_i);
            A_j                                = globalA(Z_i,Z_i);
            invATimesProjection                = A_j\diff;
            invATimesProjectionFullVector      = zeros(height,1);
            invATimesProjectionFullVector(Z_i) = invATimesProjection;
            
            
            % Old - we optimized it. Really sped-up the code.
%             Z_i      = getAdjacentEdges(i, edge, numOfNodes);
%             ZDim     = length(Z_i);
%             A_j      = zeros(ZDim, ZDim);
%             fMinusAU = globalB - globalA * U_i ;
% 
%             % Fill in A_j
%             for z = 1:ZDim
%                 for y = 1:ZDim
%                     A_j(z,y) = globalA(Z_i(z), Z_i(y));
%                 end
%             end
% 
%             fMinusAUAfterQ        = zeros(ZDim,1);
%             for j = 1:ZDim
%                 fMinusAUAfterQ(j) = fMinusAU(Z_i(j));
%             end
% 
%             invATimesProjection                       = inv(A_j) * fMinusAUAfterQ;
%             invATimesProjectionFullVector             = zeros(height,1);
%             for o = 1:ZDim
%                 invATimesProjectionFullVector(Z_i(o)) = invATimesProjection(o);
%             end
            

            U_i = U_i + invATimesProjectionFullVector;
        end

        
        % Test GS
        if(numOfGSIterations == 1)
           newval1                          = getHrCurlNormforProblem3(U_i,globalA);
           normU_iVector(numOfGSIterations) = newval1;
           errors2(numOfGSIterations)       = normU_iVector(numOfGSIterations) / norm_0;
        end
        if(numOfGSIterations > 1)
            newval2       = getHrCurlNormforProblem3(U_i,globalA);
            normU_iVector = [normU_iVector;newval2];
            newval3       = normU_iVector(numOfGSIterations) / normU_iVector(numOfGSIterations-1);
            errors2       = [errors2;newval3];
        end
        
        checkForTolerance = normU_iVector(numOfGSIterations) / norm_0;
    end    % For "Test GS".
    
    solution = U_i;
    
    % Test GS
    %To check if this is working.
    errorConvergenceRate     = 0;
    for l = 1:numOfGSIterations
        errorConvergenceRate = errorConvergenceRate + errors2(l);
    end
    errorConvergenceRate     = errorConvergenceRate/numOfGSIterations;
    
end