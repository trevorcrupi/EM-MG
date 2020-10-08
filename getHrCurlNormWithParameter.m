function error = getHrCurlNormWithParameter(c,globalA)

    error = c' * globalA * c;
    error = sqrt(abs(error));

end