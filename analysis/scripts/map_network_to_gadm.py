# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import pathlib
import datetime as dt
import pypsa
import geopandas as gpd


def parse_inputs(base_path):
    """
    The function parses the necessary inputs for the analysis
    """

    # load the data
    network_pypsa_earth_path = pathlib.Path(
        base_path, snakemake.input.pypsa_earth_network_path
    )
    gadm_shapes_path = pathlib.Path(base_path, snakemake.input.gadm_shapes_path)
    gadm_shapes = gpd.read_file(gadm_shapes_path)
    network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)

    return (
        network_pypsa_earth,
        gadm_shapes,
    )


def cluster_and_map_network(pypsa_network, gadm_dataframe):
    buses_gdf = (
        gpd.GeoDataFrame(
            pypsa_network.buses,
            geometry=gpd.points_from_xy(pypsa_network.buses.x, pypsa_network.buses.y),
            crs="EPSG:4326",
        )
        .to_crs(3857)
        .reset_index()
    )

    spatial_join_gadm_bus_gdf = (
        buses_gdf.sjoin(gadm_dataframe, how="left")
        .loc[:, ("Bus", "ISO_1")]
        .rename(columns={"ISO_1": "state_code"})
    )

    pypsa_network.generators["state"] = pypsa_network.generators.index.map(
        spatial_join_gadm_bus_gdf.to_dict()
    )

    return pypsa_network


if __name__ == "__main__":
    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "map_network_to_gadm")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "map_network_to_gadm")
    output_path = pathlib.Path(
        default_path, "analysis", "outputs", "map_network_to_gadm"
    )
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_map_to_gadm_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    (
        network_pypsa_earth_df,
        gadm_shapes_df,
    ) = parse_inputs(default_path)

    new_network = cluster_and_map_network(network_pypsa_earth_df, gadm_shapes_df)
    new_network.export_to_netcdf(
        pathlib.Path(output_path, snakemake.output.mapped_network_output_file_name)
    )

    log_output_file.close()
