function norm = getHrCurlNormForParameterParallelogram(c,iterationNum,numOfTriangles,p,t,edge,k)

    norm_u    = getL2NormParallelogram(c,iterationNum,numOfTriangles,p,t,edge,k);
    norm_curl = getL2NormCurlParallelogram(c,iterationNum,numOfTriangles,p,t,edge,k);
    
    norm = sqrt(norm_u+norm_curl);
end