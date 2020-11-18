function norm = getL2NormCurl(c,iterationNum,numOfTriangles,p,t,edge,k)

    mystr = ['newEle/new_ele' num2str(iterationNum) '.mat'];
    load(mystr);
    
    curlRZKforPhiLine1 = @(bCoeff) [ - bCoeff   ];
    curlRZKforPhiLine2 = @(bCoeff) [ - bCoeff/k ];
    curlRZKforPhiLine3 = @(aCoeff) [   aCoeff   ];           
    curlRZKforPsiLine1 = @(coeffBasis, X, Y) [ - coeffBasis(3) - coeffBasis(1).*X          ];
    curlRZKforPsiLine2 = @(coeffBasis, X, Y) [ - coeffBasis(3)/k - 3.*(coeffBasis(1)/k).*X ];
    curlRZKforPsiLine3 = @(coeffBasis, X, Y) [   coeffBasis(2) - coeffBasis(1).*Y          ];
    
    
    numOfNodes     = size(p,2);
    norm           = 0;    
    for i = 1:numOfTriangles      % Calculate (complete L2-error)^2
        
        columnVector     = t(1:3, i); % Vector with points of the triangle
        rowVector        = new_ele(i, 1:3); % Vector with edges of the triangle
        localPhiCoeffs   = getLocalPhiCoeffs(p,columnVector);
        localPsiCoeffs   = getLocalPsiCoeffs(p,rowVector,edge);
        updatedRowVector = rowVector + numOfNodes; % Now numbered for n-nodes + N-edges -- bumps up numbering of edges, so node numbers are first.
        cAndRVector      = [columnVector; updatedRowVector'];
    
        phi1 = localPhiCoeffs(:,1);
        phi2 = localPhiCoeffs(:,2);
        phi3 = localPhiCoeffs(:,3);
        psi1 = localPsiCoeffs(:,1);
        psi2 = localPsiCoeffs(:,2);
        psi3 = localPsiCoeffs(:,3);
        
        curl_rzk_line1 = @(X, Y) c(cAndRVector(1))*curlRZKforPhiLine1(phi1(2)) + c(cAndRVector(2))*curlRZKforPhiLine1(phi2(2)) + c(cAndRVector(3))*curlRZKforPhiLine1(phi3(2)) + c(cAndRVector(4))*curlRZKforPsiLine1(psi1, X, Y) + c(cAndRVector(5))*curlRZKforPsiLine1(psi2, X, Y) + c(cAndRVector(6))*curlRZKforPsiLine1(psi3, X, Y);
        curl_rzk_line2 = @(X, Y) c(cAndRVector(1))*curlRZKforPhiLine2(phi1(2)) + c(cAndRVector(2))*curlRZKforPhiLine2(phi2(2)) + c(cAndRVector(3))*curlRZKforPhiLine2(phi3(2)) + c(cAndRVector(4))*curlRZKforPsiLine2(psi1, X, Y) + c(cAndRVector(5))*curlRZKforPsiLine2(psi2, X, Y) + c(cAndRVector(6))*curlRZKforPsiLine2(psi3, X, Y);
        curl_rzk_line3 = @(X, Y) c(cAndRVector(1))*curlRZKforPhiLine3(phi1(1)) + c(cAndRVector(2))*curlRZKforPhiLine3(phi2(1)) + c(cAndRVector(3))*curlRZKforPhiLine3(phi3(1)) + c(cAndRVector(4))*curlRZKforPsiLine3(psi1, X, Y) + c(cAndRVector(5))*curlRZKforPsiLine3(psi2, X, Y) + c(cAndRVector(6))*curlRZKforPsiLine3(psi3, X, Y);
    
        [X,Y,Wx,Wy]         = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]);        
        
        localL2ErrorSquared = Wx' * (X .* (curl_rzk_line1(X,Y).^2 + curl_rzk_line2(X,Y).^2 + curl_rzk_line3(X,Y).^2)) * Wy;
        norm                = norm + localL2ErrorSquared;
    end
end