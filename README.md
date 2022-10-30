A group of classes to implement the Bayesian Optimisation algorithm. Several
surrogate models and acquisition function approaches are implemented.

Requires version 2022a of Matlab or later:
Required toolboxes: optimisation, statistics and machine learning

List of classes:

bayesOpt            - Master class. Implements the Bayesian Optimisation algorithm.
acqFcn              - An abstract class defining the acquisition function interface and
                      a template method to maximise the function.
expectedImprovement - Concrete EI acquisition function implementation
ucb                 - Concrete upper confidence bound acquisition function implementation
surrogateModel      - Abstract surrogate model interface.
gpr                 - Concrete gaussian process regression model. Wrapper for RegressionGP class
rf                  - Concrete random forest model. Wrapper for TreeBagger class
gprFitMethod        - Enumeration class for gpr model fitting methods
gprPredMethod       - Enumeration class for gpr prediction methods
kernels             - Enumeration class for supported gaussian process kernel functions
ei                  - Expected Improvement acquisition function implementation


List of live scripts

BayesOptExample.mlx - A simple example maximising the function x.sin(x) in the interval 0<=x<=10.
                      Shows how the objects and classes are combined to maximise the target
                      function. Tested 11/08/2022
