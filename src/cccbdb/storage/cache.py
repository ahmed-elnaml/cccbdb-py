"""Content-addressable local disk cache for downloaded NIST CCCBDB HTML pages."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Dict, Optional


class DiskCache:
    """File-based local cache storing raw HTML pages and provenance metadata by URL SHA-256."""

    def __init__(self, cache_dir: Path, ttl_seconds: Optional[int] = 60 * 60 * 24 * 30) -> None:
        self.cache_dir = Path(cache_dir)
        self.pages_dir = self.cache_dir / "pages"
        self.pages_dir.mkdir(parents=True, exist_ok=True)
        self.ttl_seconds = ttl_seconds

    @staticmethod
    def key_for_url(url: str) -> str:
        """Compute deterministic SHA-256 cache key for a canonical URL."""
        return hashlib.sha256(url.strip().encode("utf-8")).hexdigest()

    def _paths(self, url: str) -> tuple[Path, Path]:
        key = self.key_for_url(url)
        return self.pages_dir / f"{key}.html", self.pages_dir / f"{key}.meta.json"

    def get(self, url: str) -> Optional[str]:
        """Return cached HTML string for ``url`` if present and not expired, else ``None``."""
        html_path, meta_path = self._paths(url)
        if not html_path.exists():
            return None
        if self.ttl_seconds is not None and meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                cached_at = float(meta.get("cached_at", 0.0))
                if (time.time() - cached_at) > self.ttl_seconds:
                    return None
            except (ValueError, OSError):
                pass
        return html_path.read_text(encoding="utf-8")

    def set(self, url: str, html: str) -> Path:
        """Write HTML content and metadata for ``url`` to the local disk cache."""
        html_path, meta_path = self._paths(url)
        html_path.write_text(html, encoding="utf-8")
        meta: Dict[str, Any] = {
            "url": url,
            "cached_at": time.time(),
            "sha256": hashlib.sha256(html.encode("utf-8", errors="replace")).hexdigest(),
            "size_bytes": len(html.encode("utf-8", errors="replace")),
        }
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return html_path

    def clear(self) -> int:
        """Remove all cached HTML and metadata files. Returns number of cached entries deleted."""
        count = 0
        for p in self.pages_dir.glob("*.html"):
            p.unlink(missing_ok=True)
            count += 1
        for p in self.pages_dir.glob("*.meta.json"):
            p.unlink(missing_ok=True)
        return count

    def stats(self) -> Dict[str, Any]:
        """Return cache statistics (entry count and total size in bytes)."""
        html_files = list(self.pages_dir.glob("*.html"))
        total_bytes = sum(f.stat().st_size for f in html_files)
        return {
            "cache_dir": str(self.cache_dir),
            "entry_count": len(html_files),
            "total_bytes": total_bytes,
        }
