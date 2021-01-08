function norm = getHrCurlNormForParameterLShapeD(c,iterationNum,numOfTriangles,p,t,edge,k)

    norm_u    = getL2NormLShapeD(c,iterationNum,numOfTriangles,p,t,edge,k);
    norm_curl = getL2NormCurlLShapeD(c,iterationNum,numOfTriangles,p,t,edge,k);
    
    norm = sqrt(norm_u+norm_curl);
end