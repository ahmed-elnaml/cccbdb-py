"""Configuration management for cccbdb-py."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

__version__ = "0.1.0"
PARSER_VERSION = "0.1.0"

DEFAULT_BASE_URL = "https://cccbdb.nist.gov"
DEFAULT_USER_AGENT = (
    f"cccbdb-py/{__version__} (+https://github.com/your-username/cccbdb-py; "
    "Computational Chemistry Research Client)"
)


@dataclass
class CCCBDBConfig:
    """Runtime configuration for :class:`cccbdb.CCCBDB` client and transport layer.

    Parameters
    ----------
    base_url : str
        Root URL of the NIST CCCBDB server (default: ``https://cccbdb.nist.gov``).
    cache_enabled : bool
        Whether to cache downloaded HTML pages and normalized records locally.
    cache_dir : Path
        Directory path for local disk cache and SQLite storage.
    cache_ttl_seconds : Optional[int]
        Cache expiration time in seconds (default: 30 days; ``None`` means never expire).
    offline : bool
        If ``True``, never perform outbound network requests and raise :class:`CacheMissError`
        when a resource is not found in local cache or fixtures.
    rate_limit_delay : float
        Minimum polite delay in seconds between consecutive HTTP requests to NIST (default: 1.5s).
    timeout_seconds : float
        HTTP request timeout in seconds (default: 30.0s).
    max_retries : int
        Maximum number of automatic retries on transient HTTP 5xx or connection errors.
    user_agent : str
        Custom HTTP User-Agent header identifying the scientific client.
    strict_validation : bool
        Whether to enforce strict unit, stoichiometry, and provenance validation on parsed records.
    """

    base_url: str = field(
        default_factory=lambda: os.environ.get("CCCBDB_BASE_URL", DEFAULT_BASE_URL)
    )
    cache_enabled: bool = field(
        default_factory=lambda: os.environ.get("CCCBDB_CACHE", "1").lower() in ("1", "true", "yes")
    )
    cache_dir: Path = field(
        default_factory=lambda: Path(
            os.environ.get("CCCBDB_CACHE_DIR", str(Path.home() / ".cache" / "cccbdb-py"))
        )
    )
    cache_ttl_seconds: Optional[int] = 60 * 60 * 24 * 30  # 30 days
    offline: bool = field(
        default_factory=lambda: os.environ.get("CCCBDB_OFFLINE", "0").lower() in ("1", "true", "yes")
    )
    rate_limit_delay: float = field(
        default_factory=lambda: float(os.environ.get("CCCBDB_DELAY", "1.5"))
    )
    timeout_seconds: float = field(
        default_factory=lambda: float(os.environ.get("CCCBDB_TIMEOUT", "30.0"))
    )
    max_retries: int = 3
    user_agent: str = DEFAULT_USER_AGENT
    strict_validation: bool = True

    @property
    def sqlite_path(self) -> Path:
        """Return the path to the local normalized SQLite database."""
        return self.cache_dir / "cccbdb_normalized.sqlite"
