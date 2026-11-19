"""Project access to the supported AdaptEnv model.

Scientific analysis helpers live in their explicit submodules. Historical
calibration workflows require their recorded source revision and are no longer
advertised as the project's model API.
"""
from importlib import import_module
from ._compat import ensure_external_paths

__all__ = [
    "ThermodynamicConfig", "ThermodynamicLandscape", "LandscapeConfig",
    "FitnessLandscapeEnv", "FitnessLandscape", "PUBLIC_MODEL_SYSTEM",
    "LandscapeService", "LandscapeServiceError", "available_presets",
    "build_landscape_config_from_preset", "build_landscape_from_preset", "resolve_landscape_preset",
    "load_project_config",
]


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(name)
    ensure_external_paths()
    value = getattr(import_module("adaptenv"), name)
    globals()[name] = value
    return value
