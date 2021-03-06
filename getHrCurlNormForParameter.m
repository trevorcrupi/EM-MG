function norm = getHrCurlNormForParameter(c,iterationNum,numOfTriangles,p,t,edge,k)

    norm_u    = getL2Norm(c,iterationNum,numOfTriangles,p,t,edge,k);
    norm_curl = getL2NormCurl(c,iterationNum,numOfTriangles,p,t,edge,k);
    
    norm = sqrt(norm_u+norm_curl);
end