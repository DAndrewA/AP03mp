# equation to calculate a weighted average and error from data points with associated errors in two seperate arrays
def weightedMean(data,error):
    mError = 0
    mean = 0
    for i in range(len(error)):
        mError += 1/(error[i]**2)
        mean += data[i]/(error[i]**2)
    mean = mean/mError
    mError = mError**(-0.5)

    return mean, mError
