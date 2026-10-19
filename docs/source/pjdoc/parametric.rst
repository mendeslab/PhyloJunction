|pj| ships with a (growing) number of parametric distributions.
Distribution constructors accept scalar literals or model expressions. A constructor returns a
distribution object; ``~`` draws from it. Use ``iid`` for vectors and ``--runs`` for complete
model executions. Distribution parameters are not automatically broadcast across vectors.
On this page you will find examples of how to invoke the available distributions (using the *phylojunction* language) when building a model.

Uniform
-------

The function for assigning a uniform (``unif``) distribution to a random variable has two parameters:

    | ``min`` (real number, required): Lower bound (inclusive).
    | ``max`` (real number, required): Upper bound (inclusive).

.. code-block:: 
    :caption: **Example command.** Drawing two IID values of a uniformly distributed random variable :math:`x`, from -1.0 and 1.0.

    x ~ iid(2, unif(min=-1.0, max=1.0))

Exponential
-----------

The function for assigning a exponential (``exponential``) distribution to a random variable has two parameters:

    | ``rate`` (positive real, required): Rate (or scale) of distribution.
    | ``rate_parameterization`` (string boolean, optional): Whether the value passed in ``rate`` is the distribution's rate. Defaults to "true".

.. code-block:: 
    :caption: **Example command.** Drawing two IID values of a exponentially (scale :math:`\frac{1}{\lambda}=1.0`) distributed random variable :math:`x`.

    x ~ iid(2, exponential(rate=1.0, rate_parameterization="false"))

Gamma
-----

The function for assigning a gamma (``gamma``) distribution to a random variable has three parameters:

    | ``shape`` (positive real, required): Shape of distribution.
    | ``scale`` (positive real, required): Scale (or rate) of distribution.
    | ``rate_parameterization`` (string boolean, optional): Whether the value passed in ``scale`` is the distribution's rate. Defaults to "false".

.. code-block:: 
    :caption: **Example command.** Drawing two IID values of a gamma-distributed (shape :math:`k=10.0`, scale :math:`\theta=1.0`) random variable :math:`x`.

    x ~ iid(2, gamma(shape=10.0, scale=1.0))

Normal
------

The function for assigning a normal (``normal``) distribution to a random variable has two parameters:

    | ``mean`` (real, required): Mean (location) of distribution.
    | ``sd`` (positive real, required): Standard deviation (scale) of distribution.

.. code-block:: 
    :caption: **Example command.** Drawing two IID values of a normally (mean :math:`\mu=0.0`, standard deviation :math:`\sigma=1.0`) distributed random variable :math:`x`.

    x ~ iid(2, normal(mean=0.0, sd=1.0))

Log-normal
----------

The function for assigning a log-normal (``lognormal``) distribution to a random variable has three parameters:

    | ``meanlog`` (real, required): Mean (location) of distribution over the logarithm of the random variable.
    | ``sdlog`` (positive real, required): Standard deviation (scale) of distribution over the logarithm of the random variable.
    | ``log_space`` (string boolean, optional): Whether the value of ``meanlog`` is in log-space. Defaults to "true".

.. code-block:: 
    :caption: **Example command.** Drawing two IID values of a log-normally (mean log(:math:`\mu)=-3.25`, standard deviation :math:`\sigma=0.25`) distributed random variable :math:`x`.

    x ~ iid(2, lognormal(meanlog=-3.25, sdlog=0.25))

Categorical
-----------

The function for assigning a categorical (``categorical``) distribution to a random variable has two parameters:

    | ``cats`` (list, required): Category values.
    | ``probs`` (list, required): Nonnegative weights, one for each category; total weight must be positive.

.. code-block::
    :caption: **Example command.** Drawing two IID values of a random variable representing to two equiprobable categories,

    x ~ iid(2, categorical(cats=[0,1], probs=[0.5, 0.5]))

Beta
----

``beta(alpha, beta)`` draws a value between zero and one. Both shape parameters must be
finite and strictly positive. The mean is ``alpha / (alpha + beta)``.

.. code-block:: text

    prior := beta(alpha=2, beta=5)
    probability ~ prior

Distribution objects and IID vectors
------------------------------------

The same constructor is used for named distributions, direct sampling, and distribution-valued
arguments. For example:

.. code-block:: text

    d := normal(0, 1)
    x ~ d
    y ~ iid(3, d)
    z ~ iid(2, iid(3, normal(0, 1)))
    selected <- z[2][1]

``iid(n, distribution)`` requires a nonnegative integer length. Zero produces an empty list.
Nested IID draws remain nested lists. Indexing is one-based, accepts integer-valued expressions,
and rejects out-of-range indices. An indexed expression depends on both its vector and its
index variables. Elements do not create separate DAG nodes.

In Python, use ``distribution.draw()`` for one value and ``distribution.draw_many(count)``
for repeated draws from the same parameters. For complete model executions, use
``compile_model(source).execute(runs=10, random_seed=42)``. The resulting ``ExecutionResults``
holds values separately from the shared graph; ``results.value("x", 0)`` selects the first
execution (Python execution indices are zero-based).
