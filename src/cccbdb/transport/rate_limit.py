"""Polite rate limiter for outbound requests to NIST CCCBDB."""

from __future__ import annotations

import threading
import time


class RateLimiter:
    """Thread-safe minimum-interval rate limiter for polite scientific data access."""

    def __init__(self, delay_seconds: float = 1.5) -> None:
        self.delay_seconds = max(0.0, float(delay_seconds))
        self._lock = threading.Lock()
        self._last_request_time: float = 0.0

    def wait(self) -> float:
        """Sleep if necessary to enforce ``delay_seconds`` since the previous request.

        Returns the actual sleep duration in seconds.
        """
        if self.delay_seconds <= 0.0:
            return 0.0
        with self._lock:
            now = time.monotonic()
            elapsed = now - self._last_request_time
            sleep_for = max(0.0, self.delay_seconds - elapsed)
            if sleep_for > 0.0:
                time.sleep(sleep_for)
            self._last_request_time = time.monotonic()
            return sleep_for
