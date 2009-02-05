import vigra

def tensor2Gradient(tensor):
    """rebuild vectors in dir. of large eigenvectors with lengths sqrt(trace)"""
    return vigra.transformImage(
        vigra.tensorEigenRepresentation(tensor),
        "\l e: sqrt(e[0]+e[1])*Vector(cos(-e[2]), sin(-e[2]))", {})

def gaussianHessian(image, sigma):
    """Return three-band tensor image containing the second
    derivatives calculated with gaussianDerivativeKernels."""
    result = vigra.PythonImage(3, image.size())
    for i in range(3):
        kernel = vigra.gaussianDerivativeKernel(sigma, 2-i, sigma, i)
        result[i].copyValues(vigra.convolveImage(image, kernel))
    return result

def gm2Gradient(gm, sigma):
    """'reverse-engineer' gradient vector image from scalar image
    via Hessian from Gaussian derivative filters of given scale"""
    
    # create vectors in direction of negative curvature:
    er = vigra.tensorEigenRepresentation(
        vigra.gaussianHessian(gm, sigma))
    return vigra.transformImage(
        er, gm, "\l e,i: i*Vector(-sin(e[2]), cos(e[2]))")

def gradientTensor(img, scale):
    """Returns a 3-band image with gradient tensor summed over all
    bands.  For single-band images, this is equivalent to::
    
      vigra.vectorToTensor(
          vigra.gaussianGradientAsVector(img, scale)

    (Otherwise, it is the sum of the same for all bands.)"""
    
    bandTensors = [
        vigra.vectorToTensor(
        vigra.gaussianGradientAsVector(img.subImage(i), scale))
        for i in range(img.bands())]
    if img.bands() > 1:
        return sum(bandTensors[1:], bandTensors[0])
    return bandTensors[0]

def colorGradient(img, scale, sqrt = True):
    """Calculate Gaussian color gradient.  Calculates sum of Gaussian
    gradient tensors in all single bands, and returns pair of (gm,
    grad) where gm is the square root of the trace of the color
    gradient tensor and grad is a 2-band image with the large
    eigenvectors of the appropriate sqrt(gm2) lengths.

    If `sqrt` is set to False, returns (gm2, grad) instead, i.e. does
    not apply sqrt to the first image (slight optimization if you
    don't need the gm image)."""
    
    colorTensor = gradientTensor(img, scale)
    # gm2 := sum of squared magnitudes of grad. in each channel
    gm2 = vigra.tensorTrace(colorTensor)
    # rebuild vector in dir. of large eigenvector with length sqrt(gm2):
    grad = vigra.transformImage( # FIXME: use tensor2Gradient?
        vigra.tensorEigenRepresentation(colorTensor), gm2,
        "\l e, mag2: sqrt(mag2)*Vector(cos(-e[2]), sin(-e[2]))", {})
    if sqrt:
        gm = vigra.transformImage(gm2, "\l x: sqrt(x)")
        return gm, grad
    return gm2, grad

def gaussianGradient(img, scale, sqrt = True):
    """Return (gm, grad) tuple, which is either the result of
    `colorGradient`, or of the usual `vigra.gaussianGradient` family
    for single-band images."""
    
    if img.bands() > 1:
        return colorGradient(img, scale, sqrt)

    grad = vigra.gaussianGradientAsVector(img, scale)
    if sqrt:
        gm = vigra.norm(grad)
    else:
        gm = vigra.transformImage(grad, "\l x: squaredNorm(x)")
    return gm, grad

def structureTensor(img, innerScale, outerScale):
    """Returns structure tensor as 3-band image.  Equivalent to
    vigra.gaussianSmoothing(gradientTensor(img, innerScale),
    outerScale)."""
    return vigra.gaussianSmoothing(gradientTensor(img, innerScale), outerScale)

def structureTensorFromGradient(grad, scale):
    return vigra.gaussianSmoothing(vigra.vectorToTensor(grad), scale)
