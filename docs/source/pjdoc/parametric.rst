|pj| provides the following parametric distributions for use in model scripts.

Distribution objects and IID vectors
------------------------------------

Construct a distribution with ``:=`` and draw from it with ``~``. Use ``iid`` for repeated
independent draws with the same parameters:

.. code-block:: text

    d := normal(0, 1)
    x ~ d
    y ~ iid(3, d)
    z ~ iid(2, iid(3, d))
    selected <- z[2][1]

Here ``x`` is a scalar, ``y`` has three elements, and ``z`` contains two three-element lists.
``iid(n, distribution)`` requires a nonnegative integer length; zero produces an empty list.
Indices start at one and must be integer-valued and within bounds.
Use :ref:`model replicates <executions>` to redraw the whole model, or
:ref:`vectorization and comprehensions <vectorization>` to vary parameters between elements.

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

Cauchy
------

``cauchy(location, scale)`` describes a real-valued distribution with heavy tails.
Both parameters are required: ``location`` is a finite median and mode, and ``scale``
is finite and strictly positive. Scale is not a standard deviation; the Cauchy
has no finite mean or variance. Sample averages need not stabilize with more draws.

.. code-block:: text

    initial := cauchy(location=0, scale=1)
    x ~ initial
    xs ~ iid(10, initial)

Use it as ``start_distribution=initial`` in ``quasse`` to redraw the origin trait
for each biological attempt. ``examples/quasse_cauchy_origin.pj`` demonstrates this.
Cauchy draws also export to RevBayes as ``dnCauchy(location=..., scale=...)``.

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

``beta(alpha, beta)`` constructs a distribution on values between zero and one.
Both shape parameters must be finite and strictly positive. The mean is ``alpha / (alpha + beta)``.

.. code-block:: text

    prior := beta(alpha=2, beta=5)
    probability ~ prior
