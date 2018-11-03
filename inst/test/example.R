###### examples in Vignette ######
library(owea)
library(gtools)

## For CrossOver Dropout ##
# example 1 #
# p = 4, t = 4, n = 16, drop mechanism = (0,0,0.5,0.5)
# opt = 0, D-optimal; opt = 1, A-optimal.
set.seed(232)

# D-optimal Design
example1 <- design('dropout', n = 16, opt = 0, t = 4, p = 4,
                   drop = c(0, 0, 0.5, 0.5), max_iter = 40)

summary(example1) # printing output

eff(example1) # for efficiency

design_compare <- cbind(t(matrix(c(2,4,3,3,1,4,2,2,2,3,1,1,3,4,1,1,3,1,
                                   2,2,4,1,3,3,3,2,4,4,2,1,4,4,1,2,3,4,
                                   1,2,3,4,1,3,4,2,2,4,1,3,4,3,2,1,4,3,
                                   2,1,4,2,1,3,3,1,4,2), ncol=16)),1)


eff(example1, ex = design_compare) # relative efficiency

effLB(example1) # for lower bound efficiency

# A-optimal Design
example1 <- design('dropout', n = 16, opt = 1, t = 4, p = 4,
                   drop = c(0, 0, 0.5, 0.5), max_iter = 40)

summary(example1) # printing output

eff(example1) # for efficiency

eff(example1, ex = design_compare) # relative efficiency

effLB(example1) # for lower bound efficiency

# Example 2 #
# t = 4, p = 4, n = 19, dropout mechanism = (0,0,0.5,0.5)
# opt = 0, D-optimal; opt = 1, A-optimal.

set.seed(232)

# D-optimal Design
example2 <- design('dropout', n = 16, opt = 0, t = 4, p = 4,
                   drop = c(0, 0, 0.5, 0.5), max_iter = 40)

summary(example2) # printing output

eff(example2) # for efficiency

effLB(example2) # for lower bound efficiency

# A-optimal Design
example2 <- design('dropout', n = 16, opt = 1, t = 4, p = 4,
                   drop = c(0, 0, 0.5, 0.5), max_iter = 40)

summary(example2) # printing output

eff(example2) # for efficiency

effLB(example2) # for lower bound efficiency




## For CrossOver Proportional ##
# Example 1 #
# n = 36, t = 3, p = 3, sigma = diag(1,p),
# tau = matrix(sqrt(1+3),nrow=3, ncol=1), lambda = 0.2
# All designs are locally optimal 
set.seed(123)

# designs in literature #
design_compare <- matrix(rep(c(1,1,2,2,3,3,2,3,1,3,1,2,3,2,3,1,2,1),each=6),ncol=3)
design_compare <- cbind(design_compare,1)


# D-optimal Design #
example1 <- design('proportional', n = 36, opt = 0, t = 3, p = 3, sigma = diag(1,3),
                   tau = matrix(sqrt(1+3),nrow=3, ncol=1), lambda = 0.2,
                   max_iter = 20)

summary(example1)

eff(example1)

eff(example1, design_compare)


# A-optimal Design #
example1 <- design('proportional', n = 36, opt = 1, t = 3, p = 3, sigma = diag(1,3),
                   tau = matrix(sqrt(1+3),nrow=3, ncol=1), lambda = 0.2,
                   max_iter = 20)

summary(example1)

eff(example1)

eff(example1, design_compare)



# Example 2 #
# n = 20, t = 4, p = 3, sigma = diag(1,p),
# tau = matrix(sqrt(1+3),nrow=3, ncol=1), lambda = 0.2
# All designs are locally optimal 

# D-optimal Design #
example2 <- design('proportional', n = 20, opt = 0, t = 4, p = 3, sigma = diag(1,3),
                   tau = matrix(sqrt(1+4),nrow=4, ncol=1), lambda = 0.2,
                   max_iter = 20)

summary(example2)

eff(example2)


# A-optimal Design #
example2 <- design('proportional', n = 20, opt = 1, t = 4, p = 3, sigma = diag(1,3),
                   tau = matrix(sqrt(1+4),nrow=4, ncol=1), lambda = 0.2,
                   max_iter = 20)

summary(example2)

eff(example2)








## For Interference Model ##
# Example 1 #
# n = 10, t = 4, p = 4, sigma = diag(1,p),
set.seed(456)

# designs in literature #
design_compare <- matrix(c(2,1,4,3,1,1,3,2,4,3,2,1,4,3,2,4,4,3,2,2
                           ,1,3,3,1,4,3,2,1,1,4,1,4,2,2,4,3,2,4,3,1),ncol=4)
design_compare <- cbind(design_compare,1)


# D-optimal Design #
example1 <- design('interference', n = 10, opt = 0, t = 4, p = 4, sigma = diag(1,4),
                   max_iter = 40)

summary(example1)

eff(example1)

eff(example1, design_compare)


# A-optimal Design #
example1 <- design('interference', n = 10, opt = 1, t = 4, p = 4, sigma = diag(1,4),
                   max_iter = 20)

summary(example1)

eff(example1)

eff(example1, design_compare)



# Example 2 #
# n = 24, t = 4, p = 5, sigma = diag(1,p),

# Zheng's Universal Optimal design

design_compare <- matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,3,4,2,3,1,4,4,2,1,2,3,1,1,
                           1,1,2,2,2,3,3,3,4,4,4,2,3,4,4,3,1,2,1,4,3,1,2,4,2,
                           3,1,4,3,1,4,2,1,2,3,4,2,3,1,4,3,1,4,2,1,2,3,2,3,4,
                           4,3,1,2,1,4,3,1,2,1,1,1,2,2,2,3,3,3,4,4,4,3,4,2,3,
                           1,4,4,2,1,2,3,1,1,1,1,2,2,2,3,3,3,4,4,4),ncol=5)

design_compare <- cbind(design_compare, 1)




# D-optimal Design #
example2 <- design('interference', n = 24, opt = 0, t = 4, p = 5, sigma = diag(1,5),
                   max_iter = 50)

summary(example2)

eff(example2)

eff(example2, design_compare)

# A-optimal Design #


example2 <- design('interference', n = 24, opt = 1, t = 4, p = 5, sigma = diag(1,5),
                   max_iter = 40)

summary(example2)

eff(example2)

eff(example2, design_compare)


