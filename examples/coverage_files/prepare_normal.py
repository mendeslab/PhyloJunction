"""Write exact posterior summaries for the normal Coverage example."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from phylojunction.interface.cmdbox.cmd_parse import compile_model

example_dir = Path(__file__).resolve().parent
results = compile_model((example_dir / "normal.pj").read_text()).execute(runs=100, random_seed=42)

# For mu ~ Normal(0, 1) and y | mu ~ Normal(mu, 1), completing the square gives
# mu² + (y - mu)² = 2(mu - y/2)² + y²/2, so mu | y ~ Normal(y/2, variance=1/2).
# Normal symmetry makes the central 95% interval the HPD interval. Each row uses
# the observation from the corresponding replicate; no sorting or pooling is needed.
means = np.asarray(results.values("y")) / 2
half_width = norm.ppf(0.975) / np.sqrt(2)
pd.DataFrame({
    "posterior_mean": means,
    "lower_95hpd": means - half_width,
    "higher_95hpd": means + half_width,
}).to_csv(example_dir / "normal.csv", index=False)
print(f"Wrote {example_dir / 'normal.csv'}")
