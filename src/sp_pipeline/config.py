"""Configuration management for SP-Pipeline."""

import os
from pathlib import Path
from typing import Any, Optional

import yaml
from platformdirs import user_cache_dir, user_config_dir


APP_NAME = "sp-pipeline"
DEFAULT_CONFIG_PATH = Path(__file__).parent.parent.parent / "config" / "default_config.yaml"


def load_config(config_path: Optional[str] = None) -> dict[str, Any]:
    """Load configuration from YAML file.

    Priority:
    1. Explicit config_path argument
    2. User config dir (~/.config/sp-pipeline/config.yaml)
    3. Default config bundled with the package

    Args:
        config_path: Optional path to a custom config file.

    Returns:
        Configuration dictionary.
    """
    # Start with defaults
    config = _load_yaml(DEFAULT_CONFIG_PATH)

    # Override with user config if exists
    user_config = Path(user_config_dir(APP_NAME)) / "config.yaml"
    if user_config.exists():
        user_overrides = _load_yaml(user_config)
        config = _deep_merge(config, user_overrides)

    # Override with explicit config if provided
    if config_path:
        explicit = _load_yaml(Path(config_path))
        config = _deep_merge(config, explicit)

    # Apply environment variable overrides
    config = _apply_env_overrides(config)

    # Set computed defaults
    config = _set_defaults(config)

    return config


def _load_yaml(path: Path) -> dict:
    """Load a YAML file."""
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def _deep_merge(base: dict, override: dict) -> dict:
    """Deep merge two dictionaries. Override wins on conflicts."""
    result = base.copy()
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def _apply_env_overrides(config: dict) -> dict:
    """Apply environment variable overrides.

    Supports:
        SP_PIPELINE_NCBI_EMAIL
        SP_PIPELINE_NCBI_API_KEY
        SP_PIPELINE_SIGNALP_TOKEN  (BioLib API token for SignalP predictions)
    """
    env_email = os.environ.get("SP_PIPELINE_NCBI_EMAIL")
    if env_email:
        config.setdefault("ncbi", {})["email"] = env_email

    env_key = os.environ.get("SP_PIPELINE_NCBI_API_KEY")
    if env_key:
        config.setdefault("ncbi", {})["api_key"] = env_key

    env_signalp_token = os.environ.get("SP_PIPELINE_SIGNALP_TOKEN")
    if env_signalp_token:
        config.setdefault("predictors", {}).setdefault("signalp", {})["api_token"] = env_signalp_token

    return config


def _set_defaults(config: dict) -> dict:
    """Set computed default values."""
    # Auto-set NCBI rate limit based on API key
    ncbi = config.get("ncbi", {})
    if ncbi.get("api_key"):
        ncbi["rate_limit"] = 10

    # Set cache directory if not specified
    cache = config.get("cache", {})
    if not cache.get("directory"):
        cache["directory"] = str(Path(user_cache_dir(APP_NAME)) / "query_cache")
    config["cache"] = cache

    return config


def get_preset(config: dict, preset_name: str) -> dict:
    """Get a preset configuration by name.

    Args:
        config: Full configuration dictionary.
        preset_name: Name of the preset (e.g., 'human_type1').

    Returns:
        Preset configuration dictionary.

    Raises:
        ValueError: If preset not found.
    """
    presets = config.get("presets", {})
    if preset_name not in presets:
        available = ", ".join(presets.keys())
        raise ValueError(
            f"Unknown preset '{preset_name}'. Available presets: {available}"
        )
    return presets[preset_name]


def list_presets(config: dict) -> list[dict]:
    """List all available presets with their descriptions and required sources.

    Returns:
        List of dicts with 'name', 'description', and 'sources' keys.
    """
    presets = config.get("presets", {})
    return [
        {
            "name": name,
            "description": preset.get("description", ""),
            "sources": preset.get("sources", ["uniprot"]),
        }
        for name, preset in presets.items()
    ]
