# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


if config["enable"].get("final_adjustment",False):

    rule final_adjustment_overnight:
        params:
            chp_extendable_DE=config_provider("sector","chp_extendable_DE"),
        input:
            network=RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            ntc="data/TYNDP_NTC.csv",
        output:
            network=RESULTS
            + "prenetworks-adjusted/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        log:
            logs(RESULTS
            + "logs/final_adjustment_overnight_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/final_adjustment.py"

rule solve_sector_network:
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
            RESULTS
            + "prenetworks-adjusted/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
            if config["enable"].get("final_adjustment",False)
            else RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
    output:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
