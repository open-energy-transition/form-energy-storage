# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,


rule cluster_networks:
    input:
        expand(
            resources("networks/base_s_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_elec_networks:
    input:
        expand(
            resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
def create_kpi_path(w):
    return {
        f"{kpi_fn}": expand(RESULTS
        + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-"
        + f"{kpi_fn}"
        + "_{planning_horizons}.pdf",
        **config["scenario"],
        run=config["run"]["name"],
        )
        for kpi_fn in config_provider("kpi")(w)["custom_plots"]
    }

rule plot_KPIs_all:
    input:
        unpack(create_kpi_path)

rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/base_s_{clusters}_elec_l{ll}_{opts}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS + "figures/.statistics_plots_base_s_{clusters}_elec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_base_s_{clusters}_elec_l{ll}_{opts}",
            **config["scenario"],
            run=config["run"]["name"],
            kind=["production", "prices", "cross_border"],
        ),
