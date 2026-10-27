"""Stateful HTTP client managing ASP session cookies, retries, and caching for NIST CCCBDB."""

from __future__ import annotations

from typing import Dict, Optional
from urllib.parse import urlencode, urljoin

import httpx

from cccbdb.config import CCCBDBConfig
from cccbdb.exceptions import CacheMissError, SessionLostError, TransportError
from cccbdb.storage.cache import DiskCache
from cccbdb.transport.rate_limit import RateLimiter


class CCCBDBTransport:
    """HTTP transport layer for NIST CCCBDB with ASP session state and local caching.

    NIST CCCBDB uses classic ASP session cookies (`ASPSESSIONID...`) for pages such as
    ``energy2x.asp`` and ``vibs2x.asp``. When accessing a molecule-dependent page,
    this transport ensures the molecule context is initialized via ``exp2x.asp`` or
    ``getformx.asp`` and automatically recovers if the server responds with
    ``Sorry I lost your molecule``.
    """

    def __init__(
        self,
        config: Optional[CCCBDBConfig] = None,
        cache: Optional[DiskCache] = None,
        http_client: Optional[httpx.Client] = None,
    ) -> None:
        self.config = config or CCCBDBConfig()
        self.cache = cache or (
            DiskCache(self.config.cache_dir, ttl_seconds=self.config.cache_ttl_seconds)
            if self.config.cache_enabled
            else None
        )
        self.rate_limiter = RateLimiter(self.config.rate_limit_delay)
        self._external_client = http_client is not None
        self._client = http_client or httpx.Client(
            base_url=self.config.base_url,
            timeout=self.config.timeout_seconds,
            headers={"User-Agent": self.config.user_agent},
            follow_redirects=True,
        )
        self._active_species_key: Optional[str] = None

    def build_url(self, endpoint: str, params: Optional[Dict[str, object]] = None) -> str:
        """Construct the canonical absolute URL for a CCCBDB endpoint and query parameters."""
        base = self.config.base_url.rstrip("/") + "/" + endpoint.lstrip("/")
        if not params:
            return base
        sorted_items = [(k, str(v)) for k, v in sorted(params.items()) if v is not None]
        return f"{base}?{urlencode(sorted_items)}" if sorted_items else base

    def get_page(
        self,
        endpoint: str,
        params: Optional[Dict[str, object]] = None,
        species_context: Optional[str] = None,
        charge: int = 0,
        force_refresh: bool = False,
    ) -> tuple[str, str]:
        """Fetch an HTML page from cache or NIST CCCBDB, returning ``(canonical_url, html_text)``.

        Parameters
        ----------
        endpoint : str
            Relative ASP endpoint (e.g., ``exp2x.asp``, ``energy2x.asp``, ``energy3x.asp``).
        params : dict, optional
            URL query parameters (e.g., ``{"casno": "7732185", "charge": 0}``).
        species_context : str, optional
            If the endpoint relies on server-side ASP session state (like ``energy2x.asp``),
            pass the ``casno`` so the cache key is scoped to the molecule and session recovery works.
        charge : int
            Net molecular charge for session initialization.
        force_refresh : bool
            Bypass local cache and re-download from NIST CCCBDB.
        """
        effective_params = dict(params or {})
        cache_params = dict(effective_params)
        if species_context and "casno" not in cache_params:
            cache_params["_ctx_casno"] = species_context
            cache_params["_ctx_charge"] = charge

        canonical_url = self.build_url(endpoint, cache_params)
        request_url = self.build_url(endpoint, effective_params)

        if self.cache and not force_refresh:
            cached_html = self.cache.get(canonical_url)
            if cached_html is not None:
                return canonical_url, cached_html

        if self.config.offline:
            raise CacheMissError(
                f"Offline mode is enabled and '{canonical_url}' is not present in the local cache."
            )

        # Ensure ASP session has the active molecule when requesting session-scoped pages
        if species_context:
            target_key = f"{species_context}:{charge}"
            if self._active_species_key != target_key:
                self._prime_species_session(species_context, charge)

        html = self._request_with_retries("GET", request_url)

        # Detect classic NIST CCCBDB expired session message
        if "Sorry I lost your molecule" in html or "Molecule problem. Defaulted to" in html:
            if species_context:
                self._prime_species_session(species_context, charge)
                html = self._request_with_retries("GET", request_url)
            if "Sorry I lost your molecule" in html:
                raise SessionLostError(
                    f"NIST CCCBDB lost ASP session state while requesting '{request_url}'."
                )

        if self.cache:
            self.cache.set(canonical_url, html)

        return canonical_url, html

    def _prime_species_session(self, casno: str, charge: int = 0) -> None:
        """Initialize the server-side ASP session cookie for a specific molecule."""
        prime_url = self.build_url("exp2x.asp", {"casno": casno, "charge": charge})
        self._request_with_retries("GET", prime_url)
        self._active_species_key = f"{casno}:{charge}"

    def _request_with_retries(self, method: str, url: str) -> str:
        """Execute an HTTP request with rate limiting and exponential backoff retries."""
        last_exc: Optional[Exception] = None
        for attempt in range(1, self.config.max_retries + 1):
            self.rate_limiter.wait()
            try:
                response = self._client.request(method, url)
                response.raise_for_status()
                return response.text
            except httpx.HTTPError as exc:
                last_exc = exc
                if attempt == self.config.max_retries:
                    break
        raise TransportError(f"HTTP request to '{url}' failed after {self.config.max_retries} attempts: {last_exc}")

    def close(self) -> None:
        """Close the underlying HTTP session if owned by this transport."""
        if not self._external_client:
            self._client.close()

    def __enter__(self) -> "CCCBDBTransport":
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        self.close()
