|pj| ships with a (growing) number of parametric distributions.
All distributions support replication (i.e., plating in the DAG) and take as arguments both direct scalar values (e.g., "1.0") and values stored in DAG nodes.
On this page you will find examples of how to invoke the available distributions (using the *phylojunction* language) when building a model.

Uniform
-------

The function for assigning a uniform (``uniform``) distribution to a random variable has four parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``min`` (real number, required): Lower bound (inclusive).
    | ``max`` (real number, required): Upper bound (inclusive).

.. code-block:: 
    :caption: **Example command.** Specifying two samples of a twice replicated uniformly distributed random variable :math:`x`, from -1.0 and 1.0.

    x ~ unif(n=2, nr=2, min=-1.0, max=1.0)

Exponential
-----------

The function for assigning a exponential (``exponential``) distribution to a random variable has four parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``rate`` (positive real, required): Rate (or scale) of distribution.
    | ``rate_parameterization`` (string boolean, optional): Whether the value passed in ``rate`` is the distribution's rate. Defaults to "true".

.. code-block:: 
    :caption: **Example command.** Specifying two samples of a twice replicated exponentially (scale :math:`\frac{1}{\lambda}=1.0`) distributed random variable :math:`x`.

    x ~ exponential(n=2, nr=2, rate=1.0, rate_parameterization="false")

Gamma
-----

The function for assigning a gamma (``gamma``) distribution to a random variable has five parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``shape`` (positive real, required): Shape of distribution.
    | ``scale`` (positive real, required): Scale (or rate) of distribution.
    | ``rate_parameterization`` (string boolean, optional): Whether the value passed in ``scale`` is the distribution's rate. Defaults to "false".

.. code-block:: 
    :caption: **Example command.** Specifying two samples of a twice replicated gamma-distributed (shape :math:`k=10.0`, scale :math:`\theta=1.0`) random variable :math:`x`.

    x ~ gamma(n=2, nr=2, shape=10.0, scale=1.0)

Normal
------

The function for assigning a normal (``normal``) distribution to a random variable has four parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``mean`` (real, required): Mean (location) of distribution.
    | ``sd`` (positive real, required): Standard deviation (scale) of distribution.

.. code-block:: 
    :caption: **Example command.** Specifying two samples of a twice replicated normally (mean :math:`\mu=0.0`, standard deviation :math:`\sigma=1.0`) distributed random variable :math:`x`.

    x ~ normal(n=2, nr=2, mean=0.0, sd=1.0)

Log-normal
----------

The function for assigning a log-normal (``lognormal``) distribution to a random variable has five parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``meanlog`` (real, required): Mean (location) of distribution over the logarithm of the random variable.
    | ``sdlog`` (positive real, required): Standard deviation (scale) of distribution over the logarithm of the random variable.
    | ``log_space`` (string boolean, optional): Whether the value of ``mean`` is in log-space. Defaults to "true".

.. code-block:: 
    :caption: **Example command.** Specifying two samples of a twice replicated log-normally (mean log(:math:`\mu)=-3.25`, standard deviation :math:`\sigma=0.25`) distributed random variable :math:`x`.

    x ~ lognormal(n=1000, mean=-3.25, sd=0.25)

Categorical
-----------

The function for assigning a categorical (``categorical``) distribution to a random variable has four parameters:

    | ``n`` (integer, optional): Number of samples to draw (samples are i.i.d.). Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates to draw per sample. Defaults to 1.
    | ``cats`` (integer, required): Number of discrete categories in each draw (must be greater than 1).
    | ``probs`` (real number, required): Probability of each category (must add up to 1.0).

.. code-block::
    :caption: **Example command.** Specifying two samples of a twice replicated random variable representing to two equiprobable categories,

    x ~ categorical(n=2, nr=2, cats=[0,1], probs=[0.5, 0.5])
