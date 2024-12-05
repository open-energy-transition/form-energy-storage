# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


if config["enable"].get("final_adjustment",False):

    rule final_adjustment:
        input:
            network=resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
            ntc="data/TYNDP_NTC.csv",
        output:
            network=resources("networks-adjusted/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/final_adjustment.py"

rule solve_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=lambda w: (
            resources("networks-adjusted/base_s_{clusters}_elec_l{ll}_{opts}.nc")
            if config["enable"].get("final_adjustment",False)
            else resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc")
        ),
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_l{ll}_{opts}.nc",
        config=RESULTS + "configs/config.base_s_{clusters}_elec_l{ll}_{opts}.yaml",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_network/base_s_{clusters}_elec_l{ll}_{opts}_solver.log"
        ),
        python=RESULTS
        + "logs/solve_network/base_s_{clusters}_elec_l{ll}_{opts}_python.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network/base_s_{clusters}_elec_l{ll}_{opts}")
    threads: solver_threads
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        "shallow"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_operations_network:
    params:
        options=config_provider("solving", "options"),
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_l{ll}_{opts}.nc",
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_operations_network/base_s_{clusters}_elec_l{ll}_{opts}_op_solver.log"
        ),
        python=RESULTS
        + "logs/solve_operations_network/base_s_{clusters}_elec_l{ll}_{opts}_op_python.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_operations_network/base_s_{clusters}_elec_l{ll}_{opts}"
        )
    threads: 4
    resources:
        mem_mb=(lambda w: 10000 + 372 * int(w.clusters)),
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        "shallow"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_operations_network.py"
