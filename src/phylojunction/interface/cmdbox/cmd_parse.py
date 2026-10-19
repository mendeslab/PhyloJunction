"""Compile PhyloJunction scripts and execute complete models.

Compilation records dependencies without drawing values. Model executions are
explicit; iid describes vectors within one execution. Interactive additions use
ExecutionResults.append and publish the returned result only after success.
"""
from pathlib import Path
from phylojunction.pgm.pgm import DirectedAcyclicGraph

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# Compile expressions without sampling; callers choose the number of whole-model runs.
def compile_model(source, source_name="<script>"):
    model = DirectedAcyclicGraph()
    for number, line in enumerate(source.splitlines(), 1):
        model.add_definition(line, source_name, number)
    return model


# Read a script (or source text), compile its shared definitions, and execute it.
def script2dag(script_file_path_or_model_spec, in_pj_file=True, random_seed=None, runs=1):
    source = Path(script_file_path_or_model_spec).read_text() if in_pj_file else script_file_path_or_model_spec
    source_name = str(script_file_path_or_model_spec) if in_pj_file else "<script>"
    return compile_model(source, source_name).execute(runs, random_seed)


if __name__ == "__main__":
    # Run from the repository or an installed environment:
    # python -m phylojunction.interface.cmdbox.cmd_parse [example] [executions]
    # These scripts can also be copied into the GUI. Repeated complete models belong
    # to execute(runs=...); iid makes a vector within one model execution.
    import sys
    import tempfile
    import json

    examples = {
        'assignments': (
            'a <- 1\n'
            'b <- 2\n'
            'c <- [a,b,3]\n'
            'd <- c\n'
            'e := c[2]\n'
            'f <- [[1,2],[3,4]]\n'
            'g := f[2][1]'
        ),
        'distributions': (
            'u ~ iid(3,unif(0,1))\n'
            'n ~ normal(0,1)\n'
            'l ~ lognormal(0,1)\n'
            'e ~ exponential(1)\n'
            'g ~ gamma(2,1)\n'
            'b ~ beta(2,3)'
        ),
        # Paired parameter lists vectorize; a shared scalar sd is not a singleton list.
        'vectorized': (
            'm <- [0.0,0.01,0.011]\n'
            'a ~ lognormal(meanlog=m,sdlog=1.0)\n'
            'b ~ normal(m,1.0)'
        ),
        # Each execution redraws the hyperparameters; its ten observations share those values.
        'hierarchical': (
            'rv1 ~ lognormal(0.0,0.25)\n'
            'rv2 ~ lognormal(-1.5,0.25)\n'
            'rv3 ~ iid(10,normal(rv1,rv2),clamp=true)'
        ),
        'yule': (
            'l0rate := sse_rate(name="lambda", value=1.0, event="w_speciation")\n'
            'stash := sse_stash(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)\n'
            'tr ~ iid(2,discrete_sse(stash=stash, start_state=0, stop="age", stop_value=2.0, '
            'origin="true", cond_spn="false", cond_surv="true"))'
        ),
        'bisse': (
            'l0rate := sse_rate(name="lambda0", value=1.0, states=[0,0,0], event="w_speciation")\n'
            'l1rate := sse_rate(name="lambda1", value=1.0, states=[1,1,1], event="w_speciation")\n'
            'm0rate := sse_rate(name="mu0", value=0.8, states=[0], event="extinction")\n'
            'm1rate := sse_rate(name="mu1", value=0.8, states=[1], event="extinction")\n'
            'q01rate := sse_rate(name="q01", value=0.6, states=[0,1], event="transition")\n'
            'q10rate := sse_rate(name="q10", value=0.6, states=[1,0], event="transition")\n'
            'stash := sse_stash(flat_rate_mat=[l0rate, l1rate, m0rate, m1rate, q01rate, q10rate], '
            'n_states=2, n_epochs=1)\n'
            'tr ~ discrete_sse(stash=stash, start_state=0, stop="age", stop_value=1.0, origin="true", '
            'cond_spn="false", cond_surv="true")'
        ),
        # Time-heterogeneous Yule process: four epochs ending at ages 2.2, 1.2, and 0.7.
        'time_yule': (
            'l0ratet1 := sse_rate(name="lambdat1", value=1.0, event="w_speciation", epoch=1)\n'
            'l0ratet2 := sse_rate(name="lambdat2", value=0.25, event="w_speciation", epoch=2)\n'
            'l0ratet3 := sse_rate(name="lambdat3", value=3.0, event="w_speciation", epoch=3)\n'
            'l0ratet4 := sse_rate(name="lambdat4", value=0.4, event="w_speciation", epoch=4)\n'
            'stash := sse_stash(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3, l0ratet4], n_states=1, '
            'n_epochs=4, seed_age=3.0, epoch_age_ends=[2.2, 1.2, 0.7])\n'
            'tr ~ iid(2,discrete_sse(stash=stash, start_state=0, stop="age", stop_value=3.0, '
            'origin="true", cond_spn="false", cond_surv="true"))'
        ),
        'fossils': (
            'birth_rate := sse_rate(name="lambda", value=1.0, states=[0,0,0], event="w_speciation")\n'
            'death_rate := sse_rate(name="mu", value=0.5, states=[0], event="extinction")\n'
            'fossil_rate := sse_rate(name="psi", value=0.8, states=[0], event="anc_sampling")\n'
            'stash := sse_stash(flat_rate_mat=[birth_rate,death_rate,fossil_rate], n_states=1, '
            'n_epochs=1)\n'
            'trs ~ iid(5,discrete_sse(stash=stash, start_state=0, stop="age", stop_value=3.0, '
            'origin="true", cond_spn="true", cond_surv="true"))'
        ),
        # Time-heterogeneous birth/death rates, with the boundary at age 1.5.
        'time_birth_death': (
            'birth_rate_t0 <- 0.5\n'
            'death_rate_t0 <- 0.25\n'
            'birth_rate_t1 <- 1.0\n'
            'death_rate_t1 <- 0.4\n'
            'det_birth_rate_t0 := sse_rate(name="lambda_t0", value=birth_rate_t0, states=[0,0,0], '
            'event="w_speciation", epoch=1)\n'
            'det_death_rate_t0 := sse_rate(name="mu_t0", value=death_rate_t0, states=[0], '
            'event="extinction", epoch=1)\n'
            'det_birth_rate_t1 := sse_rate(name="lambda_t1", value=birth_rate_t1, states=[0,0,0], '
            'event="w_speciation", epoch=2)\n'
            'det_death_rate_t1 := sse_rate(name="mu_t1", value=death_rate_t1, states=[0], '
            'event="extinction", epoch=2)\n'
            'stash := sse_stash(flat_rate_mat=[det_birth_rate_t0, det_death_rate_t0, det_birth_rate_t1,'
            ' det_death_rate_t1], n_states=1, n_epochs=2, epoch_age_ends=[1.5], seed_age=3.0)\n'
            'tr ~ discrete_sse(stash=stash,start_state=0,stop="age",stop_value=3.0,origin=true)'
        ),
        # Time-heterogeneous GeoSSE: two regions plus their combined state, in two epochs.
        'time_geosse': (
            'w_birth_rate0_t0 <- 1.5\n'
            'death_rate0_t0 <- 0.25\n'
            'trans_rate02_t0 <- 0.0\n'
            'w_birth_rate1_t0 <- 1.0\n'
            'death_rate1_t0 <- 0.5\n'
            'trans_rate12_t0 <- 0.4\n'
            'b_birth_rate201_t0 <- 2.0\n'
            'b_birth_rate202_t0 <- 1.0\n'
            'b_birth_rate212_t0 <- 1.0\n'
            'trans_rate20_t0 <- 0.6\n'
            'trans_rate21_t0 <- 0.6\n'
            'w_birth_rate0_t1 <- 0.85\n'
            'death_rate0_t1 <- 0.5\n'
            'trans_rate02_t1 <- 0.8\n'
            'w_birth_rate1_t1 <- 1.5\n'
            'death_rate1_t1 <- 0.75\n'
            'trans_rate12_t1 <- 0.6\n'
            'b_birth_rate201_t1 <- 2.5\n'
            'b_birth_rate202_t1 <- 1.25\n'
            'b_birth_rate212_t1 <- 1.25\n'
            'trans_rate20_t1 <- 0.75\n'
            'trans_rate21_t1 <- 0.75\n'
            'det_w_birth_rate0_t0 := sse_rate(name="lambda0_t0", value=w_birth_rate0_t0, '
            'states=[0,0,0], event="w_speciation", epoch=1)\n'
            'det_death_rate0_t0 := sse_rate(name="mu0_t0", value=death_rate0_t0, states=[0], '
            'event="extinction", epoch=1)\n'
            'det_trans_rate02_t0 := sse_rate(name="q02_t0", value=trans_rate02_t0, states=[0,2], '
            'event="transition", epoch=1)\n'
            'det_w_birth_rate1_t0 := sse_rate(name="lambda1_t0", value=w_birth_rate1_t0, '
            'states=[1,1,1], event="w_speciation", epoch=1)\n'
            'det_death_rate1_t0 := sse_rate(name="mu1_t0", value=death_rate1_t0, states=[1], '
            'event="extinction", epoch=1)\n'
            'det_trans_rate12_t0 := sse_rate(name="q12_t0", value=trans_rate12_t0, states=[1,2], '
            'event="transition", epoch=1)\n'
            'det_b_birth_rate201_t0 := sse_rate(name="lambda201_t0", value=b_birth_rate201_t0, '
            'states=[2,0,1], event="bw_speciation", epoch=1)\n'
            'det_b_birth_rate202_t0 := sse_rate(name="lambda202_t0", value=b_birth_rate202_t0, '
            'states=[2,0,2], event="bw_speciation", epoch=1)\n'
            'det_b_birth_rate212_t0 := sse_rate(name="lambda212_t0", value=b_birth_rate212_t0, '
            'states=[2,1,2], event="bw_speciation", epoch=1)\n'
            'det_trans_rate20_t0 := sse_rate(name="q20_t0", value=trans_rate20_t0, states=[2,0], '
            'event="transition", epoch=1)\n'
            'det_trans_rate21_t0 := sse_rate(name="q21_t0", value=trans_rate21_t0, states=[2,1], '
            'event="transition", epoch=1)\n'
            'det_w_birth_rate0_t1 := sse_rate(name="lambda0_t1", value=w_birth_rate0_t1, '
            'states=[0,0,0], event="w_speciation", epoch=2)\n'
            'det_death_rate0_t1 := sse_rate(name="mu0_t1", value=death_rate0_t1, states=[0], '
            'event="extinction", epoch=2)\n'
            'det_trans_rate02_t1 := sse_rate(name="q02_t1", value=trans_rate02_t1, states=[0,2], '
            'event="transition", epoch=2)\n'
            'det_w_birth_rate1_t1 := sse_rate(name="lambda1_t1", value=w_birth_rate1_t1, '
            'states=[1,1,1], event="w_speciation", epoch=2)\n'
            'det_death_rate1_t1 := sse_rate(name="mu1_t1", value=death_rate1_t1, states=[1], '
            'event="extinction", epoch=2)\n'
            'det_trans_rate12_t1 := sse_rate(name="q12_t1", value=trans_rate12_t1, states=[1,2], '
            'event="transition", epoch=2)\n'
            'det_b_birth_rate201_t1 := sse_rate(name="lambda201_t0", value=b_birth_rate201_t1, '
            'states=[2,0,1], event="bw_speciation", epoch=2)\n'
            'det_b_birth_rate202_t1 := sse_rate(name="lambda202_t0", value=b_birth_rate202_t1, '
            'states=[2,0,2], event="bw_speciation", epoch=2)\n'
            'det_b_birth_rate212_t1 := sse_rate(name="lambda212_t0", value=b_birth_rate212_t1, '
            'states=[2,1,2], event="bw_speciation", epoch=2)\n'
            'det_trans_rate20_t1 := sse_rate(name="q20_t1", value=trans_rate20_t1, states=[2,0], '
            'event="transition", epoch=2)\n'
            'det_trans_rate21_t1 := sse_rate(name="q21_t1", value=trans_rate21_t1, states=[2,1], '
            'event="transition", epoch=2)\n'
            'stash := sse_stash(flat_rate_mat=[det_w_birth_rate0_t0, det_death_rate0_t0, '
            'det_trans_rate02_t0, det_w_birth_rate1_t0, det_death_rate1_t0, det_trans_rate12_t0, '
            'det_b_birth_rate201_t0, det_b_birth_rate202_t0, det_b_birth_rate212_t0, '
            'det_trans_rate20_t0, det_trans_rate21_t0, det_w_birth_rate0_t1, det_death_rate0_t1, '
            'det_trans_rate02_t1, det_w_birth_rate1_t1, det_death_rate1_t1, det_trans_rate12_t1, '
            'det_b_birth_rate201_t1, det_b_birth_rate202_t1, det_b_birth_rate212_t1, '
            'det_trans_rate20_t1, det_trans_rate21_t1], n_states=3, n_epochs=2, epoch_age_ends=[1.0], '
            'seed_age=3.0)\n'
            'trs ~ iid(2,discrete_sse(stash=stash, start_state=0, stop="age", stop_value=3.0, '
            'origin="true"))'
        ),
        # Construct a stash for each rate explicitly; each stash produces two IID trees.
        'incomplete_sampling': (
            'birth_rate <- [0.8,0.9]\n'
            'det_birth_rate := '
            'sse_rate(name="lambda",value=birth_rate,states=[0,0,0],event="w_speciation")\n'
            'det_sampling_rate := sse_prob(name="rho",value=0.5,state=0)\n'
            'stashes := [sse_stash(flat_rate_mat=[rate],flat_prob_mat=[det_sampling_rate],n_states=1) '
            'for rate in det_birth_rate]\n'
            'trs ~ [iid(2,discrete_sse(stash=s,stop="age",stop_value=1.0,origin=true)) for s in '
            'stashes]'
        ),
        # read_tree always returns a list; [1] selects its first annotated tree.
        'read_tree': (
            'tr <- '
            'read_tree(string="((sp1[&index=1]:1.0,sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)[&index=5];",'
            ' node_name_attr="index")\n'
            'first := tr[1]'
        ),
        'read_tree_unannotated': (
            'tr <- read_tree(string="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")\n'
            'first := tr[1]'
        ),
    }
    # File reading uses two embedded trees, so this example also works outside the checkout.
    examples['read_tree_file'] = ''
    # Invalid examples are deliberately separate: missing parameters and unequal vector
    # lengths should explain the problem instead of producing partial samples.
    errors = {
        'error_missing': 'u ~ unif(max=1)',
        'error_lengths': 'a ~ lognormal(meanlog=[0,0.01],sdlog=[1,0.9,0.8])',
    }
    example = sys.argv[1] if len(sys.argv) > 1 else 'hierarchical'
    if example not in examples and example not in errors:
        raise SystemExit('Choose an example: ' + ', '.join([*examples, *errors]))
    runs = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    source = examples[example] if example in examples else errors[example]
    if example == 'read_tree_file':
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'trees.tre'
            path.write_text('((sp1:1,sp2:1):1,sp3:2);\n((sp4:1,sp5:1):1,sp6:2);\n')
            source = f'trees <- read_tree(file_path={json.dumps(str(path))})\nsecond := trees[2]'
            results = compile_model(source, '<example:read_tree_file>').execute(runs, random_seed=42)
    else:
        results = compile_model(source, f'<example:{example}>').execute(runs, random_seed=42)
    for name in results.model.name_node_dict:
        print(name, results.values(name))
