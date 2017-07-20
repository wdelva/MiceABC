dummy.wrapper <- function(inputvector = input.vector){ # Let's start with 4 input parameters, e.g. 100, 2, 3, 4
  random.gamma.data <- rgamma(n = inputvector[1], shape = inputvector[2], scale = inputvector[3])
  random.norm.data <- rnorm(n = inputvector[1], mean = inputvector[2]*inputvector[3], sd = inputvector[4])
  feature1 <- median(random.gamma.data)
  feature2 <- sd(random.gamma.data * random.norm.data)
  feature3 <- IQR(random.gamma.data + random.norm.data)
  feature4 <- diff(range(random.gamma.data * random.norm.data))

  outputvector <- c(feature1, feature2, feature3, feature4)

  return(outputvector)
}
