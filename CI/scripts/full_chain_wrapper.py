import alice3_full_chain
import job_configs.ChainConfig

import argparse
import os
import sys
import pathlib
import importlib


def load_config(path: str):

    config_dir = str(pathlib.Path(path).parent.resolve())
    if config_dir not in sys.path:
        sys.path.insert(0, config_dir)

    spec = importlib.util.spec_from_file_location("job_config", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    return mod.cfg


def print_cfg(cfg):

    import dataclasses
    from pprint import pprint

    for name, value in vars(type(cfg)).items():
        if name.startswith("_"):
            continue
        if dataclasses.is_dataclass(value):
            print(f"--- {name} ---")
            pprint(dataclasses.asdict(value))
            print()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run Alice3 full chain test with optional geometry argument.")
    parser.add_argument(
        "--geometry",
        "-g",
        default="",
        help="Optional geometry argument to pass to the full chain."
    )

    args, _ = parser.parse_known_args()
    full_chain_args, _ = alice3_full_chain.getArgumentParser().parse_known_args()
    cfg = load_config(full_chain_args.config)

    if (args.geometry == ""):
        print("No geometry argument provided. Using default geometry.")
    else:
        print(f"Using geometry: {args.geometry}")
        cfg.general.geo_dir = "geometries/" + args.geometry + "/geom"
        cfg.general.digi_file = "geometries/" + args.geometry + "/digiConfigurations/digi.json"

        # check if the geometry directory exists
        if not os.path.exists(cfg.general.geo_dir):
            raise FileNotFoundError(f"Geometry directory {cfg.general.geo_dir} does not exist."
                                    f" Please check the geometry argument.")

        # Select Gen1/Gen3: Gen1 geometries ship a python builder
        # (buildALICE3Geometry.py); Gen3 ones are built by the compiled blueprint
        # builder (buildALICE3GeometryGen3) and have no such file.
        cfg.general.geomGeneration = (
            "gen1" if os.path.exists(cfg.general.geo_dir + "/buildALICE3Geometry.py")
            else "gen3")
        print(f"Geometry generation: {cfg.general.geomGeneration}")

        # if digi.json does not exist, check for any other json file and use that instead
        if not os.path.exists(cfg.general.digi_file):
            print(f"Warning: {cfg.general.digi_file} does not exist " +
                  f"Checking for other digi<>.json files in geometries/{args.geometry}/digiConfigurations...")
            list_of_files = os.listdir("geometries/" + args.geometry + "/digiConfigurations")
            json_digi_files = sorted(
                [f for f in list_of_files if f.endswith(".json") and f.startswith("digi")])
            n_files_found = len(json_digi_files)
            if n_files_found >= 1:  # one or more found: pick deterministically
                if n_files_found > 1:
                    print(
                        f"Warning: multiple digi<>.json files in geometries/{args.geometry}/digiConfigurations: "
                        f"{json_digi_files}. Using '{json_digi_files[0]}'.")
                else:
                    print(
                        f"Found {json_digi_files[0]} in geometries/{args.geometry}/digiConfigurations. Using that instead.")
                cfg.general.digi_file = "geometries/" + args.geometry + "/digiConfigurations/" + json_digi_files[0]
            else:
                raise FileNotFoundError(
                    f"No digi.json or other digi<>.json file found in geometries/{args.geometry}/digiConfigurations")

    # execute

    alice3_full_chain.runFullChain(cfg, args=full_chain_args)
