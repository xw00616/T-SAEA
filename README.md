# T-SAEA
Despite the success of evolutionary algorithms (EAs) for solving multi-objective problems, most of them are based on the
assumption that all objectives can be evaluated within the
same period of time. However, in many real-world applications, such an assumption is unrealistic since different objectives
must be evaluated using different computer simulations or
physical experiments with various time complexities. To address this issue, a surrogate assisted evolutionary algorithm
along with a parameter-based transfer learning (T-SAEA)
is proposed in this work. While the surrogate for the cheap
objective can be updated on sufficient training data, the surrogate for the expensive one is updated by either the training
data set or a transfer learning approach. To find out the
transferable knowledge, a filter-based feature selection algorithm is used to capture the pivotal features of each objective,
and then use the common important features as a carrier for
knowledge transfer between the cheap and expensive objectives. Then, the corresponding parameters in the surrogate
models are adaptively shared to enhance the quality of the
surrogate models.

# Reference
If you found T-SAEA useful, we would be grateful if you cite the following reference:

Wang, Xilu, Yaochu Jin, Sebastian Schmitt, and Markus Olhofer. "Transfer learning for gaussian process assisted evolutionary bi-objective optimization for objectives with different evaluation times." In Proceedings of the 2020 Genetic and Evolutionary Computation Conference, pp. 587-594. 2020.
